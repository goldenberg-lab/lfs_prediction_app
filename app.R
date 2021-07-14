library(shiny)
library(shinymaterial)

# overall help here: https://ericrayanderson.github.io/shinymaterial/
# color universe here: http://materializecss.com/color.html
# icon universe here :http://materializecss.com/icons.html

ui <- material_page(
  title = "LFS prediction",
  nav_bar_color = "blue",
  nav_bar_fixed = FALSE,
  
  
  # Place side-nav in the beginning of the UI
  material_side_nav(fixed = FALSE, # Place side-nav tabs within side-nav
                    material_side_nav_tabs(side_nav_tabs = c("Example tab 1" = "example_1",
                                                             "Example tab 2" = "example_2",
                                                             "About" = "about_section"),
                                           icons = c("assessment", "blur_off","favorite_border" ))
  ),
  
  material_row(
    material_column(
      width = 2,
      material_card(
        title = "",
        depth = 4,
        material_switch(
          input_id = "example_input",
          off_label = "",
          on_label = ""
        ))),
    material_column(width = 10,
                    # Define side-nav tab content
                    material_side_nav_tab_content(
                      side_nav_tab_id = "example_1",
                      tags$h4("First tab")
                    ),
                    
                    # Define side-nav tab content
                    material_side_nav_tab_content(
                      side_nav_tab_id = "example_2",
                      tags$h4("Second tab")
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
  
}
shinyApp(ui = ui, server = server)