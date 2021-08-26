library(shiny)
library(shinymaterial)
library(data.table)

# overall help here: https://ericrayanderson.github.io/shinymaterial/
# color universe here: http://materializecss.com/color.html
# icon universe here :http://materializecss.com/icons.html

source('global.R')
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
      width = 2,
      material_card(
        title = "Card title",
        depth = 4,
        fileInput('target_upload', 'Choose file to upload',
                  accept = c(
                    'text/csv',
                    'text/comma-separated-values',
                    '.csv',
                    '.rda'
                  ))
        )),
    material_column(width = 10,
                    # Define side-nav tab content
                    material_side_nav_tab_content(
                      side_nav_tab_id = "pred",
                      DT::dataTableOutput("sample_table"),
                      plotOutput('author_plot_null'),
                      plotOutput('author_plot_cancer')
                
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
  df_products_upload <- reactive({
    inFile <- input$target_upload
    if (is.null(inFile))
      return(NULL)
    df <- read.csv(inFile$datapath)
    # df <- df[, 1:20]
    # save(df, file = 'test.rda')
    # df <- read.csv(inFile$datapath, nrows = 2)
    df$gender <- ifelse(df$gender == "M", 0, 1)
    df$cancer_atdraw <- ifelse(df$cancer_atdraw == 'Yes', 1, 0)
    df$systemic_treatment_atdraw <- ifelse(df$systemic_treatment_atdraw == 'Yes', 1, 0)
    test_y <- factor(ifelse(df$ageofonset > 72 | is.na(df$ageofonset), "No", "Yes"))
    clin_names <- names(df)[!names(df) %in% feat_list]
    clin_dat <- df[, clin_names]
    df_mat <- df[,!names(df) %in% clin_names]
    df_mat <- scale(df_mat)
    test.predictions <- predict(mod
                                ,newdata = df_mat
                                ,type='prob')
    
    valid_results <- as.data.frame(cbind(test_pred = valid.predictions, test_label = valid_y, valid_clin))
    
    test_results <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))
    
    # run model
    
    
    return(df)
  })
  
  output$sample_table<- DT::renderDataTable({
    df <- df_products_upload()
    DT::datatable(df)
  })
  
  output$author_plot_null <- renderPlot({
    # plot
    get_results(temp = test, age_pred = 6, show_metrics = FALSE)[[2]]
 
  })
  
  output$author_plot_cancer <- renderPlot({
    # plot
  
    get_results(temp = test, age_pred = 6, show_metrics = FALSE)[[1]]
  })
}
shinyApp(ui = ui, server = server)