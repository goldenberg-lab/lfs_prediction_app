library(shiny)
library(shinymaterial)
library(data.table)

# overall help here: https://ericrayanderson.github.io/shinymaterial/
# color universe here: http://materializecss.com/color.html
# icon universe here :http://materializecss.com/icons.html

# increase limit for data upload to 30MB
options(shiny.maxRequestSize = 30*1024^2)

source('global.R')
source('utils.R')
ui <- material_page(
  title = "LFS prediction",
  nav_bar_color = "blue",
  nav_bar_fixed = FALSE,
  
  # Place side-nav in the beginning of the UI
  material_side_nav(fixed = FALSE, # Place side-nav tabs within side-nav
                    material_side_nav_tabs(side_nav_tabs = c("Predictions" = "pred",
                                                             "About" = "about_section"),
                                           icons = c("assessment", "blur_off" ))
  ),
  
  material_row(
    material_column(
      width = 3,
      material_card(
        title = "1. Upload your data",
        depth = 4,
        fileInput('target_upload', 'Choose file to upload',
                  accept = c(
                    'text/csv',
                    'text/comma-separated-values',
                    '.csv',
                    '.rda'
                  ))
        ),
      br(), br(),
      
      material_card(
        title = '2. Apply data corrections',
        depth = 4,
        material_checkbox(input_id = 'remove_outliers',
                          label = 'Remove outliers',
                          initial_value = FALSE),
        br(), br(),
        material_checkbox(input_id = 'array_correction',
                          label = 'PC correction for array (only if your data has 450k)',
                          initial_value = FALSE),
        br(), br(),
        material_checkbox(input_id = 'batch_correction',
                          label = 'Batch correction?',
                          initial_value = FALSE)
      )
      ),
    material_column(width = 9,
                    # Define side-nav tab content
                    material_side_nav_tab_content(
                      side_nav_tab_id = "pred",
                      material_card(
                        depth = 4,
                        uiOutput('directions_and_title'),
                        DT::dataTableOutput("sample_table"),
                        material_checkbox(input_id = 'compare_results',
                                          label = "Visualize predictions with the author's results",
                                          initial_value = FALSE),
                      ),
                      material_column(width = 12,
                                      plotOutput('author_plot')
                                      
                      )
                     
                    ),
                    
                    # Define side-nav tab content
                    material_side_nav_tab_content(
                      side_nav_tab_id = "about_section",
                      tags$h4("About tab"),
                      
                      # reference webiste 
                      tags$a(
                        target = "_blank",
                        class = "btn orange",
                        href = "goldenberglab.ca",
                        "Visite our website!"
                      )
                    )
    )
  )    
)

server <- function(input, output) {
  
  # reactive data frame that the user uploads. 
  user_preds <- reactive({
    ac <- input$array_correction
    bc <- input$batch_correction
    rl <- input$remove_outliers
    
    inFile <- input$target_upload
    if (is.null(inFile))
      return(NULL)
    # This can be replaced later in the pipeline so that this part starts once removing confounders is done
    # df <- readRDS('df_small.rda')
    df <- readRDS(inFile$datapath)

    if(rl){
      message('---- Removing outliers')
      
      # remove outliers
      all_probes <- colnames(df)[grepl('cg',colnames(df))]
      all_clin <- names(df)[!names(df) %in% all_probes]
      pc <- prcomp(as.matrix(df[,all_probes]), scale = TRUE)
      pc_clin_before <- cbind(df[,all_clin],pc$x)
      keep <- get_outliers(pc_clin_before,3)
      df <- df[df$SentrixID %in% keep,]
    }
    if(ac){
      message('---- Begin array correction')
      ## Remove array confounder ## 
      df_450k <- df[df$array == "450",]
     
      ## correct 450k data ##
      corrected_450k <- remove_array_confounder(df_450k)
      df_850k <- df[df$array == "850",]
      df <- rbind(corrected_450k,df_850k)
    
    }
    
    if(bc){
      message('---- Begin batch correction')
      df <- remove_batch_confounder(df)
    }
    
    ## Read in features ## 
    features <- read.csv('features.txt',sep='\t')
    features$probe <- as.character(features$probe)
    features$gene <- as.character(features$gene)
    ## Aggregate probes  ## 
    df <- aggregate_probes(df,features)
    ## Scale data ## 
    df <- scale_df(df,features$gene)
    
    ## Predict on new data ## 
    xgboost_results <- pred_cancer_xgboost_test(df,features$gene)
    
    ## Calibrate results ## 
    calibrated_results <- platt_scaling(xgboost_results)
    ## Generate prediction metrics ## 
    ROCobj_test <- ROCInfo_atcutoff(calibrated_results,other_title = 'Predictions')
    
    # run model
    result_dat <- calibrated_results
    # save(result_dat, file ='temp_results.rda')
    
    return(result_dat)
  })
  
  output$directions_and_title <- renderUI({
    inFile <- input$target_upload
    if(is.null(inFile)){
      h5("Once you've uploaded data, the model will run automatically and show results below")
    } else {
      h5('Predicted probabilities with clinical data')
    }
  })
  
  output$sample_table<- DT::renderDataTable({
    df <- user_preds()
    if(is.null(df)){
      NULL
    } else {
      df$test_pred.No <- NULL
      names(df)[1] <- 'Predicted probability'
      names(df)[2] <- 'Ground truth' 
      DT::datatable(df)
    }
    
  })
  
  output$author_plot <- renderPlot({
    cr <- input$compare_results 
    if(cr){
      df <- user_preds()
      if(is.null(df)){
        NULL
      } else {
        int_names <- intersect(names(df), names(test))
        df <- df[, int_names]
        test <- test[, int_names]
        df$user <- 'User'
        test$user <- 'Author'
        df <- rbind(df, test)
      }
      # plot
      compare_results(temp = df)
    } else {
      NULL
    }
  })
}
shinyApp(ui = ui, server = server)