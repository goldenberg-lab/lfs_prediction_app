get_results <- function(temp, age_pred,show_metrics=FALSE){
  temp$preds <- temp$test_pred_calibrated.Yes
  temp <- temp[!duplicated(temp$tm_donor),]
  temp$test_label <-  factor(temp$test_label, levels = c('Yes', 'No'))
  # temp$test_label_number <- ifelse(temp$test_label == 'Yes', 1, 0)
  # temp <- calibrate_probs(test_results = temp)
  # temp$preds <- probability.calibration(as.numeric(temp$test_label_number), as.numeric(round(temp$test_pred.Yes, 3)))
  temp$pred_label <- ifelse(temp$preds >= .5, 'Yes', 'No')
  temp$pred_label <- factor(temp$pred_label, levels=c('Yes', 'No'))
  
  # HERE get sensitivity and specificity
  model_info <-caret::confusionMatrix(temp$pred_label, temp$test_label)
  f1_meas <- round(caret::F_meas( temp$pred_label, temp$test_label),2)
  mod_sens <- round(model_info$byClass[[1]], 2)
  mod_spec <-round(model_info$byClass[[2]], 2)
  auc_value <- round(pROC::auc(temp$test_label, temp$preds), 2)
  temp_null <- temp[temp$cancer_diagnosis=='Unaffected',]
  temp_cancer <- temp[temp$cancer_diagnosis!='Unaffected',]
  if(show_metrics){
    null_plot = conmat_paper(data = temp_null, predict = 'preds',actual = 'test_label', cutoff = 0.5,get_plot = TRUE,other_title = paste0('auc= ', auc_value, '   sensitivity= ',mod_sens,'   specificity= ',mod_spec, ' f1=  ', f1_meas ), data_type = 'null', age_pred = age_pred)
    cancer_plot = conmat_paper(data = temp_cancer, predict = 'preds',actual = 'test_label', cutoff = 0.5,get_plot = TRUE,other_title = paste0('auc= ', auc_value, '   sensitivity= ',mod_sens,'   specificity= ',mod_spec, ' f1=  ', f1_meas ), data_type = 'cancer', age_pred = age_pred)
    
  } else {
    null_plot = conmat_paper(data = temp_null, predict = 'preds',actual = 'test_label', cutoff = 0.5,get_plot = TRUE,other_title = '', data_type = 'null', age_pred = age_pred)
    cancer_plot = conmat_paper(data = temp_cancer, predict = 'preds',actual = 'test_label', cutoff = 0.5,get_plot = TRUE,other_title = '', data_type = 'cancer',age_pred = age_pred)
    
  }
  
  return(list(null_plot, cancer_plot))
}

calibrate_probs_valli <- function(train_results, test_results) {
  actual <- "test_label"
  predict <- "test_pred.Yes"
  ll_df <- data.frame(x=train_results[[predict]],y=as.factor(train_results[[actual]]))
  model_log <-glm(y~x,data = ll_df,family = binomial)
  #predicting on the cross validation after platt scaling
  result_platt<-predict(model_log,ll_df["x"],type = "response")
  train_results$test_pred_calibrated.Yes <- result_platt
  
  # Predicting on the test dataset 
  ll_df_test<-data.frame(x=test_results[[predict]])
  result_test_platt<-predict(model_log,ll_df_test,type="response")
  test_results$test_pred_calibrated.Yes <- result_test_platt
  
  return(list(train_results,test_results))
}

conmat_paper <- function( data, predict, actual, cutoff, get_plot, other_title, data_type, age_pred)
{	
  predict <- data[[predict]]
  actual  <- relevel( as.factor( data[[actual]] ), "Yes" )
  if(data_type == 'null') {
    age <- data$agesamplecollection
    cancer_name <- data$cancer_diagnosis
    gender <- data$gender
    cancer_atdraw <- data$cancer_atdraw
    # gender <- ifelse(data$`F`==1, 'Female', 'Male')
    # cancer_atdraw <- ifelse(data$`Y`==1, 'cancer', 'no_cancer')
    result <- data.table( actual = actual, predict = predict, age = age, cancer_name = cancer_name,
                          gender=gender, cancer_atdraw=cancer_atdraw)
    if(age_pred == 6){
      result$actual <- 'No cancer before 6'
    } else if (age_pred == 5){
      result$actual <- 'No cancer before 5'
    } else {
      result$actual <- 'No cancer before 4'
    }
    
    # caculating each pred falls into which category for the confusion matrix
    result[ , type := ifelse( predict >= cutoff, 'FP', 'TN' )
            %>% as.factor() ]
    
    library(ggrepel)
    # jittering : can spread the points along the x axis 
    result <- as.data.frame(result)
    result$age <- round(result$age/12, 2)
    plot <- ggplot( result, aes( actual, predict, color = type ) ) + 
      geom_violin( fill = "grey", color = NA ) +
      geom_jitter(size = 3, show.legend = TRUE) +
      scale_color_manual(name = '',
                         values = c('red', 'blue'),
                         breaks = c( "FP", "TN"))+
      geom_hline( yintercept = cutoff, color = 'black', alpha = 0.6, linetype = 2 ) + 
      # geom_text(aes_string(label = text_name),alpha = 0.7, fontface = "bold",position=position_jitter(width = 0.49, height = 0), show.legend = FALSE)+
      # geom_vline(xintercept = 1.5, linetype = 2) +
      scale_y_continuous( limits = c(0, 1.01 ) ) + 
      guides( col = guide_legend( nrow = 2 ) ) + # adjust the legend to have two rows  
      labs(subtitle= other_title, y='Predictions' , x = '') +
      theme(text = element_text(size=10))+
      theme(plot.subtitle = element_text(size=8)) +
      ggthemes::theme_base()
    
    if(get_plot) {
      return(plot)
    } else {
      return(as.data.frame(result))
    }
    
  } else {
    age <- data$ageofonset
    
    age_diff <- data$ageofonset - data$agesamplecollection
    cancer_atdraw <- data$cancer_atdraw
    age_sample_collect <- data$agesamplecollection
    cancer_name <- data$cancer_diagnosis
    gender <- data$gender
    result <- data.table( actual = actual, predict = predict, age = age,age_diff= age_diff, cancer_name = cancer_name, gender=gender, cancer_atdraw=cancer_atdraw)
    
    
    result$actual <- ifelse(result$actual=='Yes', paste0('Cancer before ', age_pred),paste0('Cancer after ', age_pred))
    result$actual <- factor(result$actual, levels = c( paste0('Cancer before ', age_pred), paste0('Cancer after ', age_pred)))
    
    # caculating each pred falls into which category for the confusion matrix
    result[ , type := ifelse( predict >= cutoff & actual ==  paste0('Cancer before ', age_pred), "TP",
                              ifelse( predict >= cutoff & actual == paste0('Cancer after ', age_pred), "FP", 
                                      ifelse( predict <  cutoff & actual ==paste0('Cancer before ', age_pred), "FN", "TN" ) ) ) %>% as.factor() ]
    
    library(ggrepel)
    # jittering : can spread the points along the x axis 
    result <- as.data.frame(result)
    result$age <- round(result$age/12, 2)
    plot <- ggplot( result, aes( actual, predict, color = type ) ) + 
      geom_violin( fill = "grey", color = NA ) +
      geom_jitter(size = 3, show.legend = TRUE) +
      scale_color_manual(name = '',
                         values = c('darkgreen', 'darkorange','red', 'blue'),
                         breaks = c( "TP", "FN", "FP", "TN" ))+
      geom_hline( yintercept = cutoff, color = 'black', alpha = 0.6, linetype = 2 ) + 
      # geom_text(aes_string(label = text_name),alpha = 0.7, fontface = "bold",position=position_jitter(width = 0.49, height = 0), show.legend = FALSE)+
      geom_vline(xintercept = 1.5, linetype = 2) +
      scale_y_continuous( limits = c( 0, 1.01 ) ) + 
      guides( col = guide_legend( nrow = 2 ) ) + # adjust the legend to have two rows  
      labs(subtitle= other_title, y= 'Predictions', x = '' ) +
      theme(text = element_text(size=10)) +
      theme(plot.subtitle = element_text(size=8)) +
      ggthemes::theme_base()
    
    
    plot
    if(get_plot) {
      return(plot)
    } else {
      return(as.data.frame(result))
    }
    
  }
  
}
