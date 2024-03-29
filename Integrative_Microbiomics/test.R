library(shiny)
library(shinydashboard)
library(shinyjs)

ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            useShinyjs()
        ),
        mainPanel(
            box(id = "myBox", title = "Tree Output", width = '800px',
                    selectInput(inputId = "myInput", label = "my input", choices = c(letters))
                    ),
            actionButton(inputId = "button", label = "show / hide")
        )
    )
)

server <- function(input, output){

    ## observe the button being pressed
    observeEvent(input$button, {

        if(input$button %% 2 == 1){
            shinyjs::hide(id = "myBox")
        }else{
            shinyjs::show(id = "myBox")
        }
    })
}

