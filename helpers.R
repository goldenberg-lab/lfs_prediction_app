#####################################################################
############################# Libraries #############################
#####################################################################

require(reshape2)
require(minfi)
require(dplyr)
require(doParallel)
# require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
require(stringr)
require(ggplot2)
require(ROCR)
require(data.table)
require(caret)
require(gridExtra)
require(grid)
require(MLmetrics)

#####################################################################
############################# Functions #############################
#####################################################################


remove_batch_confounder <- function(data) {
  p_train <- readRDS('rds/batchEffectRemoval.rds')
  cols <- grepl( "cg" , names( data ))
  beta <- data[ , cols ]
  clin <- data[ , -cols ]
  test_pred <-predict(p_train,beta)
  Xhat_pred <- test_pred[, !(colnames(test_pred) %in% c(maxpc))] %*% t(p_train$rotation[, !(colnames(p_train$rotation) %in% c(maxpc))])
  beta_adj_pred <- scale(Xhat_pred, center = -(colMeans(beta)), scale = T)
  data_adj_test <- cbind(clin,beta_adj_pred)
  return(data_adj_test)
}

##########
# Function to batch 450k data onto the 850k space
##########

remove_array_confounder <- function(data) {
  cat("[  450k data onto 850k space ]","\n")
  models <- readRDS('rds/arrayEffectRemoval.rds')
  probes <- grepl( "cg" , names( data )); lm_all <- list()
  for (probe in probes) {
    probe_model=models[[probe]]
    new_data <- data.frame(probe_450 = data[,probe],
                           age_sample_collection = data$agesamplecollection,
                           gender = data$gender)
    lm_probe <- predict(probe_model,new_data)
    if (i %% 50000 == 0) {
      print(i)
    }
  }
  lm_all[[probe]] <-lm_probe
  corrected_meth <- do.call(cbind,lm_all)
  corrected_data <- cbind(data[,-cols],corrected_meth)
  return(corrected_data)
}

##########
# Function that identifies ids of outliers
##########

get_outliers <- function(pc,n) {
  pc$PC1 <- scale(pc$PC1)
  pc$PC2 <- scale(pc$PC2)
  u_thres_pc1 <- mean(pc$PC1) + n*sd(pc$PC1)
  l_thres_pc1 <- mean(pc$PC1) - n*sd(pc$PC1)
  u_thres_pc2 <- mean(pc$PC2) + n*sd(pc$PC2)
  l_thres_pc2 <- mean(pc$PC2) - n*sd(pc$PC2)
  pc2 <- pc[pc$PC1 < u_thres_pc1 & pc$PC1 > l_thres_pc1, ]
  pc2 <- pc2[pc2$PC2 < u_thres_pc2 & pc2$PC2 > l_thres_pc2, ]
  outlier_id <- as.character(pc$SentrixID[!pc$SentrixID %in% pc2$SentrixID])
  cat(paste0(length(outlier_id)," outliers removed","\n"))
  keep_id <- as.character(pc2$SentrixID)
  return(keep_id)
}

##########
# Function that quantifies association between PCs and confounders
##########

generate_pcsummary <- function(pc_clin) {
  pcs <- colnames(pc_clin)[grepl( "PC" , names(pc_clin) )] 
  allp_array <- list() 
  allp_batch <- list() 
  allp_agesamplecollection<- list()
  for (i in pcs) {
    allp_array[[i]] <- wilcox.test(pc_clin[,i] ~ as.factor(pc_clin$array))$p.value
    allp_batch[[i]] <- summary(aov(pc_clin[,i] ~ as.factor(pc_clin$batch)))[[1]][1,"Pr(>F)"]
    allp_batch[[i]] <- cor.test(pc_clin[,i] ~ as.factor(pc_clin$agesamplecollection))$p.value
  }
  
  allp <- data.frame(PC=pcs,
                     batch_p=as.numeric(do.call("cbind",allp_batch)),
                     array_p=as.numeric(do.call("cbind",allp_array)),
                     agesamplecollection_p=as.numeric(do.call("cbind",allp_agesamplecollection)))
  allp <- cbind(allp,data.frame(t(summary(pc)$importance)))
  allp$batch_padj <- allp$batch_p/allp$Proportion.of.Variance
  allp$array_padj <- allp$array_p/allp$Proportion.of.Variance
  allp$agesamplecollection_padj <- allp$agesamplecollection_p/allp$Proportion.of.Variance
  return(allp)
}


##########
# Function to scale methylation variables
##########

e_df <- function(data,genes) {
  int_feat <- intersect(names(data), genes)
  tmp <- data[,int_feat]
  tmp <- scale(tmp, center=TRUE,scale=TRUE)
  clin_vars <- names(data)[!names(data) %in% int_feat]
  data_scaled <- cbind(data[,clin_vars],tmp)
  return(data_scaled)
}

##########
# Function to extract probes based on location and aggregate by gene
##########

aggregate_probes <- function(data,features) {
  data <- data[,!grepl('ch', names(data))]
  allprobes <- colnames(data)[grepl('cg',colnames(data))]
  clincols <- colnames(data)[!names(data) %in% allprobes]
  feat_data <- data[c(clincols,features$probe)]
  meth <- data.frame(t(feat_data[features$probe])) 
  colnames(meth) <- as.character(feat_data$SentrixID)
  meth$gene <- features$gene[match(rownames(meth),features$probe)]
  meth <- meth %>% group_by(gene)  %>% summarise(across(everything(), list(mean)))  %>% as.data.frame()
  meth <- meth[complete.cases(meth), ]
  cat(paste0("[ Genes ] : ", dim(meth)[1],"\n"))
  rownames(meth) <- meth$gene
  feat_data <- cbind(data[clincols],t(meth[2:(length(meth))]))
  return(feat_data)
}

##########
# Function to predict cancer before a given age of onset cutoff given a xgboost model
##########

pred_cancer_xgboost_test <- function(test_dat, features) {
  
  ## Read in predictive model ##
  model <- readRDS('NoobCorrected_after_covs_beta_ProjPC2Adj_lfs_5UTR_scaled_canceratdraw_systreat_xgboost_ageofonset_model.rds')
  
  age_cutoff = 72
  # get intersection of features and real data
  model_features <- c('gender', 'cancer_atdraw','systemic_treatment_atdraw', features)
  model_features <- intersect(model_features, colnames(test_dat))
  test_dat$gender <- ifelse(test_dat$gender == "M", 0, 1)
  test_dat$cancer_atdraw <- ifelse(test_dat$cancer_atdraw == 'Yes', 1, 0)
  test_dat$systemic_treatment_atdraw <- ifelse(test_dat$systemic_treatment_atdraw == 'Yes', 1, 0)
  
  # get y
  test_y <- factor(ifelse(test_dat$ageofonset > age_cutoff | is.na(test_dat$ageofonset), "No", "Yes"))
  
  # get clinical data
  test_clin <- test_dat[1:42]
  
  # adding this. get features from model
  mod_feats <- model[[11]]$xNames
  
  # get model data
  test_dat <- test_dat[, mod_feats]
  
  test.predictions <- predict(model,
                              data.matrix(test_dat),
                              type = 'prob')
  
  # combine predictions and real labels
  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))
  return(temp_dat)
  
}

##########
# Function to calibrate probability scores
##########

platt_scaling <- function(test_results) {
  ## Read in platt scaling recalibration model ## 
  recal_model <- readRDS('NoobCorrected_after_covs_beta_ProjPC2Adj_lfs_5UTR_scaled_canceratdraw_systreat_xgboost_recalibration_model.rds')
  # predicting on the test dataset using Platt Scaling
  ll_df_test<-data.frame(x=test_results[["test_pred.Yes"]])
  result_test_platt<-predict(recal_model,ll_df_test,type="response")
  test_results$test_pred_calibrated.Yes <- result_test_platt
  return(test_results)
}

##########
# Function to calculate F1 score
##########

get_f1 <- function (data) {
  precision <- Precision(data$test_label, data$predicted_label, positive = 1)
  recall <- Recall(data$test_label, data$predicted_label, positive = 1)
  f1 <- 2*precision*recall/(precision+recall)
  return(f1)
}

##########
# Function to get ROC, F1, sensitivity and specificity at optimal cutoff in test data
##########

ROCInfo_atcutoff <- function(data,other_title) {
  
  # decision boundary 
  cutoff <- 0.5
  
  if (any(unique(data$test_label) %in% c("Yes","No"))) {
    data$test_label <- factor(ifelse(data$test_label == "Yes", 1, 0))
  }
  
  # calculate the roc values
  pred <- prediction( data$test_pred_calibrated.Yes, data$test_label )
  perf <- performance( pred, "tpr", "fpr" )
  roc_dt <- data.frame( fpr = perf@x.values[[1]], tpr = perf@y.values[[1]] )
  auc <- performance( pred, "auc" )@y.values[[1]]
  
  ## get sens/spec/f1
  data$predicted_label <- factor(ifelse(data$test_pred_calibrated.Yes >= cutoff,1,0))
  sensitivity <- sensitivity(data$predicted_label,data$test_label,positive=1)
  specificity <- specificity(data$predicted_label,data$test_label,positive=1)
  f1 <- get_f1(data)
  
  options(scipen = '999')
  
  # the main title for the plot
  sub_title <- sprintf(other_title,  "Cutoff at %.2f , AUC = %.3f", 
                       cutoff, auc )
  # plot roc curve
  roc_plot <- ggplot( roc_dt, aes( fpr, tpr ) ) + 
    geom_line( color = rgb( 0, 0, 1, alpha = 0.3 ) ) +
    geom_segment( aes( x = 0, y = 0, xend = 1, yend = 1 ), alpha = 0.8, color = "royalblue" ) + 
    labs( title = sub_title, x = "False Postive Rate", y = "True Positive Rate" ) +
    geom_hline( yintercept = sensitivity, alpha = 0.8, linetype = "dashed", color = "steelblue4" ) +
    geom_vline( xintercept = 1-specificity, alpha = 0.8, linetype = "dashed", color = "steelblue4" ) +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5))
  
  return( list( data = data,
                plot    = roc_plot, 
                pred      = pred,
                perf      = perf,
                cutoff    = cutoff, 
                auc         = auc,
                sensitivity = sensitivity, 
                specificity = specificity,
                f1 = f1) )
}

##########
# Function that plots the top 2 PCs coloured by confounders
##########

generate_pcplots <- function(pc_clin) {
  cols <- colorRampPalette(brewer.pal(n = 9, 'Set1'))(length(unique(pc$gender)))
  pdf(paste0(outdir,"Plots/",output,"_PCA_gender.pdf"),width=9,height=7)
  cancerstatusplot <- ggplot(pc_clin,aes(PC1,PC2,color=cancerstatus)) + 
    geom_point(size = 3, alpha = 0.7) +
    xlab('PC1') + ylab('PC2') +
    scale_color_manual(name = '', values = cols) + 
    geom_hline(yintercept= 0, linetype="dashed",color = "grey", size=1) +
    geom_vline(xintercept=0, linetype="dashed", color = "grey", size=1) +
    theme_minimal()
  print(cancerstatusplot)
  suppressMessages(dev.off())
  
  pdf(paste0(outdir,"Plots/",output,"_PCA_agesamplecollection.pdf"),width=9,height=7)
  ageplot <- ggplot(pc_clin,aes(PC1,PC2,color=agesamplecollection)) + 
    geom_point(size = 3, alpha = 0.7) +
    xlab('PC1') + ylab('PC2') +
    scale_color_manual(name = '', values = cols) + 
    geom_hline(yintercept= 0, linetype="dashed",color = "grey", size=1) +
    geom_vline(xintercept=0, linetype="dashed", color = "grey", size=1) +
    theme_minimal()
  print(ageplot)
  suppressMessages(dev.off())
  
  cols <- colorRampPalette(brewer.pal(n = 9, 'Set1'))(length(unique(pc$batch)))
  pc_clin$array <- factor(pc_clin$array)
  pdf(paste0(outdir,"Plots/",output,"_PCA_confounders.pdf"),width=9,height=7)
  confounderplot <- ggplot(pc_clin,aes(PC1,PC2,color=batch,shape=array)) + 
    geom_point(size = 3, alpha = 0.7) +
    xlab('PC1') + ylab('PC2') +
    scale_color_manual(name = '', values = cols) + 
    geom_hline(yintercept= 0, linetype="dashed",color = "grey", size=1) +
    geom_vline(xintercept=0, linetype="dashed", color = "grey", size=1) +
    theme_minimal()
  print(confounderplot)
  suppressMessages(dev.off())
  
}

##########
# function that calibrates the probabilites (named after valli!)
##########
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

##########
# function that visualizes user results in the same plot as author's results
##########

compare_results <- function(temp){
  temp$preds <- temp$test_pred_calibrated.Yes
  temp <- temp[!duplicated(temp$tm_donor),]
  temp$test_label <-  factor(temp$test_label, levels = c('Yes', 'No'))
  # temp$test_label_number <- ifelse(temp$test_label == 'Yes', 1, 0)
  # temp <- calibrate_probs(test_results = temp)
  # temp$preds <- probability.calibration(as.numeric(temp$test_label_number), as.numeric(round(temp$test_pred.Yes, 3)))
  temp$pred_label <- ifelse(temp$preds >= .5, 'Yes', 'No')
  temp$pred_label <- factor(temp$pred_label, levels=c('Yes', 'No'))
  
  # sensitivity and specificity
  model_info <-caret::confusionMatrix(temp$pred_label, temp$test_label)
  f1_meas <- round(caret::F_meas( temp$pred_label, temp$test_label),2)
  mod_sens <- round(model_info$byClass[[1]], 2)
  mod_spec <-round(model_info$byClass[[2]], 2)
  auc_value <- round(pROC::auc(temp$test_label, temp$preds), 2)
  temp_null <- temp[temp$cancer_diagnosis=='Unaffected',]
  temp_cancer <- temp[temp$cancer_diagnosis!='Unaffected',]
  temp_null <- evaluate_results(temp_null, type ='null')
  temp_cancer <- evaluate_results(temp_cancer, type = 'cancer')
  dat <- rbind(temp_null, temp_cancer)
  
  ggplot( dat, aes( actual, predict, color = type ) ) + 
    geom_violin( fill = "grey", color = NA ) +
    geom_jitter(size = 3, show.legend = TRUE) +
    scale_color_manual(name = '',
                       values = c('darkgreen', 'darkorange','red', 'blue', 'gold'),
                       breaks = c( "TP", "FN", "FP", "TN", 'User' ))+
    geom_hline( yintercept = 0.5, color = 'black', alpha = 0.6, linetype = 2 ) + 
    # geom_text(aes_string(label = text_name),alpha = 0.7, fontface = "bold",position=position_jitter(width = 0.49, height = 0), show.legend = FALSE)+
    geom_vline(xintercept = 1.5, linetype = 2) +
    scale_y_continuous( limits = c( 0, 1.01 ) ) + 
    guides( col = guide_legend( nrow = 2 ) ) + # adjust the legend to have two rows  
    labs(subtitle= 'Results', y= 'Predictions', x = '' ) +
    theme(text = element_text(size=10)) +
    theme(plot.subtitle = element_text(size=8)) +
    ggthemes::theme_base()
}

##########

##########
evaluate_results <- function(data, type){
  # NULL
  predict <- data[['preds']]
  actual  <- relevel( as.factor( data[['test_label']] ), "Yes" )
  if(type == 'null'){
    age <- data$agesamplecollection
    cancer_name <- data$cancer_diagnosis
    gender <- data$gender
    cancer_atdraw <- data$cancer_atdraw
    user <- data$user
    # gender <- ifelse(data$`F`==1, 'Female', 'Male')
    # cancer_atdraw <- ifelse(data$`Y`==1, 'cancer', 'no_cancer')
    result <- data.table( actual = actual, predict = predict, age = age, cancer_name = cancer_name,
                          gender=gender, cancer_atdraw=cancer_atdraw, user=user)
    result$actual <- 'No cancer before 6'
    
    # caculating each pred falls into which category for the confusion matrix
    result[ , type := ifelse( predict >= 0.5, 'FP', 'TN' )
            %>% as.character() ]
    
    result$type <- ifelse(result$user=='User', 'User', result$type)
    
    library(ggrepel)
    # jittering : can spread the points along the x axis 
    result <- as.data.frame(result)
    result$age <- round(result$age/12, 2)
  } else {
    age <- data$ageofonset
    
    cancer_atdraw <- data$cancer_atdraw
    age_sample_collect <- data$agesamplecollection
    cancer_name <- data$cancer_diagnosis
    gender <- data$gender
    user <- data$user
    
    result <- data.table( actual = actual, predict = predict, age = age, cancer_name = cancer_name, gender=gender, cancer_atdraw=cancer_atdraw, user= user)
    
    
    result$actual <- ifelse(result$actual=='Yes', paste0('Cancer before ', '6'),paste0('No cancer before ', '6'))
    result$actual <- factor(result$actual, levels = c( paste0('Cancer before ', '6'), paste0('No cancer before ', '6')))
    
    # caculating each pred falls into which category for the confusion matrix
    result[ , type := ifelse( predict >= .5 & actual ==  paste0('Cancer before ', '6'), "TP",
                              ifelse( predict >= .5 & actual == paste0('No cancer before ', '6'), "FP", 
                                      ifelse( predict <  .5 & actual ==paste0('Cancer before ', '6'), "FN", "TN" ) ) ) %>% as.factor() ]
    
    # jittering : can spread the points along the x axis 
    result <- as.data.frame(result)
    result$age <- round(result$age/12, 2)    
  }
  
  return(result)
}

# ##########
# ##########
# get_results <- function(temp, age_pred){
#   temp$preds <- temp$test_pred_calibrated.Yes
#   temp <- temp[!duplicated(temp$tm_donor),]
#   temp$test_label <-  factor(temp$test_label, levels = c('Yes', 'No'))
#   # temp$test_label_number <- ifelse(temp$test_label == 'Yes', 1, 0)
#   # temp <- calibrate_probs(test_results = temp)
#   # temp$preds <- probability.calibration(as.numeric(temp$test_label_number), as.numeric(round(temp$test_pred.Yes, 3)))
#   temp$pred_label <- ifelse(temp$preds >= .5, 'Yes', 'No')
#   temp$pred_label <- factor(temp$pred_label, levels=c('Yes', 'No'))
#   
#   # HERE get sensitivity and specificity
#   model_info <-caret::confusionMatrix(temp$pred_label, temp$test_label)
#   f1_meas <- round(caret::F_meas( temp$pred_label, temp$test_label),2)
#   mod_sens <- round(model_info$byClass[[1]], 2)
#   mod_spec <-round(model_info$byClass[[2]], 2)
#   auc_value <- round(pROC::auc(temp$test_label, temp$preds), 2)
#   temp_null <- temp[temp$cancer_diagnosis=='Unaffected',]
#   temp_cancer <- temp[temp$cancer_diagnosis!='Unaffected',]
#   
#     null_plot = conmat_paper(data = temp_null, predict = 'preds',actual = 'test_label', cutoff = 0.5,get_plot = TRUE,other_title = '', data_type = 'null', age_pred = age_pred)
#     cancer_plot = conmat_paper(data = temp_cancer, predict = 'preds',actual = 'test_label', cutoff = 0.5,get_plot = TRUE,other_title = '', data_type = 'cancer',age_pred = age_pred)
#     
#   
#   
#   return(list(null_plot, cancer_plot))
# }



##########
# function that creates the confusion matrix plot
##########
# 
# conmat_paper <- function( data, predict, actual, cutoff, get_plot, other_title, data_type, age_pred)
# {	
#   predict <- data[[predict]]
#   actual  <- relevel( as.factor( data[[actual]] ), "Yes" )
#   if(data_type == 'null') {
#     age <- data$agesamplecollection
#     cancer_name <- data$cancer_diagnosis
#     gender <- data$gender
#     cancer_atdraw <- data$cancer_atdraw
#     user <- data$user
#     # gender <- ifelse(data$`F`==1, 'Female', 'Male')
#     # cancer_atdraw <- ifelse(data$`Y`==1, 'cancer', 'no_cancer')
#     result <- data.table( actual = actual, predict = predict, age = age, cancer_name = cancer_name,
#                           gender=gender, cancer_atdraw=cancer_atdraw, user=user)
#     if(age_pred == 6){
#       result$actual <- 'No cancer before 6'
#     } else if (age_pred == 5){
#       result$actual <- 'No cancer before 5'
#     } else {
#       result$actual <- 'No cancer before 4'
#     }
#     
#     # caculating each pred falls into which category for the confusion matrix
#     result[ , type := ifelse( predict >= cutoff, 'FP', 'TN' )
#             %>% as.character() ]
#     
#     result$type <- ifelse(result$user=='User', 'User', result$type)
#     
#     library(ggrepel)
#     # jittering : can spread the points along the x axis 
#     result <- as.data.frame(result)
#     result$age <- round(result$age/12, 2)
#     plot <- ggplot( result, aes( actual, predict, color = type ) ) + 
#       geom_violin( fill = "grey", color = NA ) +
#       geom_jitter(size = 3, show.legend = TRUE) +
#       scale_color_manual(name = '',
#                          values = c('red', 'blue', 'green'),
#                          breaks = c( "FP", "TN",'User'))+
#       geom_hline( yintercept = cutoff, color = 'black', alpha = 0.6, linetype = 2 ) + 
#       # geom_text(aes_string(label = text_name),alpha = 0.7, fontface = "bold",position=position_jitter(width = 0.49, height = 0), show.legend = FALSE)+
#       # geom_vline(xintercept = 1.5, linetype = 2) +
#       scale_y_continuous( limits = c(0, 1.01 ) ) + 
#       guides( col = guide_legend( nrow = 2 ) ) + # adjust the legend to have two rows  
#       labs(subtitle= other_title, y='Predictions' , x = '') +
#       theme(text = element_text(size=10))+
#       theme(plot.subtitle = element_text(size=8)) +
#       ggthemes::theme_base()
#     
#     if(get_plot) {
#       return(plot)
#     } else {
#       return(as.data.frame(result))
#     }
#     
#   } else {
#     age <- data$ageofonset
#     
#     age_diff <- data$ageofonset - data$agesamplecollection
#     cancer_atdraw <- data$cancer_atdraw
#     age_sample_collect <- data$agesamplecollection
#     cancer_name <- data$cancer_diagnosis
#     gender <- data$gender
#     result <- data.table( actual = actual, predict = predict, age = age,age_diff= age_diff, cancer_name = cancer_name, gender=gender, cancer_atdraw=cancer_atdraw)
#     
#     
#     result$actual <- ifelse(result$actual=='Yes', paste0('Cancer before ', age_pred),paste0('Cancer after ', age_pred))
#     result$actual <- factor(result$actual, levels = c( paste0('Cancer before ', age_pred), paste0('Cancer after ', age_pred)))
#     
#     # caculating each pred falls into which category for the confusion matrix
#     result[ , type := ifelse( predict >= cutoff & actual ==  paste0('Cancer before ', age_pred), "TP",
#                               ifelse( predict >= cutoff & actual == paste0('Cancer after ', age_pred), "FP", 
#                                       ifelse( predict <  cutoff & actual ==paste0('Cancer before ', age_pred), "FN", "TN" ) ) ) %>% as.factor() ]
#     
#     library(ggrepel)
#     # jittering : can spread the points along the x axis 
#     result <- as.data.frame(result)
#     result$age <- round(result$age/12, 2)
#     plot <- ggplot( result, aes( actual, predict, color = type ) ) + 
#       geom_violin( fill = "grey", color = NA ) +
#       geom_jitter(size = 3, show.legend = TRUE) +
#       scale_color_manual(name = '',
#                          values = c('darkgreen', 'darkorange','red', 'blue'),
#                          breaks = c( "TP", "FN", "FP", "TN" ))+
#       geom_hline( yintercept = cutoff, color = 'black', alpha = 0.6, linetype = 2 ) + 
#       # geom_text(aes_string(label = text_name),alpha = 0.7, fontface = "bold",position=position_jitter(width = 0.49, height = 0), show.legend = FALSE)+
#       geom_vline(xintercept = 1.5, linetype = 2) +
#       scale_y_continuous( limits = c( 0, 1.01 ) ) + 
#       guides( col = guide_legend( nrow = 2 ) ) + # adjust the legend to have two rows  
#       labs(subtitle= other_title, y= 'Predictions', x = '' ) +
#       theme(text = element_text(size=10)) +
#       theme(plot.subtitle = element_text(size=8)) +
#       ggthemes::theme_base()
#     
#     
#     plot
#     if(get_plot) {
#       return(plot)
#     } else {
#       return(as.data.frame(result))
#     }
#     
#   }
#   
# }


