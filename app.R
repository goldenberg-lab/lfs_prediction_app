library(shiny)
library(shinymaterial)
library(data.table)

# Using shinymaterial package. links below have help.
# overall help here: https://ericrayanderson.github.io/shinymaterial/
# color universe here: http://materializecss.com/color.html
# icon universe here :http://materializecss.com/icons.html

# increase limit for data upload to 30MB
options(shiny.maxRequestSize = 30*1024^2)

source('global.R')
#source('utils.R') # commented out because now using file called helpers.R which is sourced in global.R
# helpers.R is has all the stuff that was in utils.R, I just added more functions.

# define the UI. the material_page function is the same thing as the fluidPage function in normal shiny. All the UI contents go into it. material_page is only differnt in aesthetic. 
ui <- material_page(
  title = "LFS prediction",
  nav_bar_color = "blue",
  nav_bar_fixed = FALSE,
  
  # defines the sections in the app if you click on the 3 bars on top left. Two sections: Predictions and About.
  material_side_nav(fixed = FALSE, # Place side-nav tabs within side-nav
                    material_side_nav_tabs(side_nav_tabs = c("Predictions" = "pred",
                                                             "About" = "about_section"),
                                           icons = c("assessment", "blur_off" ))
  ),
  
  # same as fluidRow in normal shiny
  material_row(
    # same as the column function in normal shiny
    material_column(
      width = 3,
      # this displays the html card around the input
      material_card(
        title = "1. Upload your data",
        depth = 4,
        # input for uploading data. First argument is the name of the object to be called on the server side (target_upload), the second is the label that shows on the UI (Choose a file to upload)
        fileInput('target_upload', 'Choose file to upload',
                  accept = c(
                    'text/csv',
                    'text/comma-separated-values',
                    '.csv',
                    '.rda'
                  ))
        ),
      # HTML breaks that provide space between inputs
      br(), br(),
      # the html card that holds the following inputs (purely aesthetic): remove_outliers, array_correction, batch_correction.
      material_card(
        title = '2. Apply data corrections',
        depth = 4,
        # input checkbox for user to remove outliers
        material_checkbox(input_id = 'remove_outliers',
                          label = 'Remove outliers',
                          initial_value = FALSE),
        br(), br(),
        # input check box for user to implement arrya correction (not working yet because need RDS folder. line 40 in global.R)
        material_checkbox(input_id = 'array_correction',
                          label = 'PC correction for array (only if your data has 450k)',
                          initial_value = FALSE),
        br(), br(),
        # input check box for user to implement batch correction (not working yet. line 25 in global.R)
        material_checkbox(input_id = 'batch_correction',
                          label = 'Batch correction?',
                          initial_value = FALSE)
      )
      ),
    # defines the column that holds the title, table and plot (right hand side of page)
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
                    # placeholder for the about section
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

# the server function takes the input ids from the UI function and performs the analysis.
server <- function(input, output) {
  
  # the reactive function is initiated when the data is chosen to upload
  user_preds <- reactive({
    # if user clicks on the checkbox for the below inputs (called using input$) then those inputs (ac, bc, rl) will be boolean objects (true, false).
    ac <- input$array_correction
    bc <- input$batch_correction
    rl <- input$remove_outliers
    
    inFile <- input$target_upload
    if (is.null(inFile))
      return(NULL)
    
    # reads in data
    df <- readRDS(inFile$datapath)

    # if the input$remove_outliers input is true, peform outlier removal
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
    # if input$array_correction is true (user click on the checkbox), then perform array correction (not working yet)
    if(ac){
      message('---- Begin array correction')
      ## Remove array confounder ## 
      df_450k <- df[df$array == "450",]
     
      ## correct 450k data ##
      corrected_450k <- remove_array_confounder(df_450k)
      df_850k <- df[df$array == "850",]
      df <- rbind(corrected_450k,df_850k)
    
    }
    
    # if input$batch_correction true, then perofrm batch correction (not working yet)
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

    return(result_dat)
  })
  
  # this output object displays text. If no file is uploaded, it displays the first h5 element. if data uploaded, it changes.
  output$directions_and_title <- renderUI({
    inFile <- input$target_upload
    if(is.null(inFile)){
      h5("Once you've uploaded data, the model will run automatically and show results below")
    } else {
      h5('Predicted probabilities with clinical data')
    }
  })
  
  # this output object is table. Its non existant initially (null), once the reactive data frame user_preds() is initiated (data uploaded), the table is rendered
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
  
  # this output object is the plot. it takes to input: (1) the input$compare_results (true of false based on if user clicked the checkbox) and (2) the reactive data frame that was uploaded with or without the corrections.
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