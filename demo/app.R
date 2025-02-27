library(shiny)
library(AIPW)
library(readr)

ui <- fluidPage(
  titlePanel("AIPW Estimation of ATE"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload CSV File", accept = ".csv"),
      uiOutput("var_select"),
      actionButton("run", "Estimate ATE")
    ),
    mainPanel(
      verbatimTextOutput("ate_output")
    )
  )
)

server <- function(input, output, session) {
  data <- reactive({
    req(input$file)
    read_csv(input$file$datapath)
  })

  output$var_select <- renderUI({
    req(data())
    selectInput("outcome", "Select Outcome Variable", choices = names(data()))
    selectInput("treatment", "Select Treatment Variable", choices = names(data()))
  })

  estimate_ate <- eventReactive(input$run, {
    req(input$outcome, input$treatment)
    df <- data()
    covariates <- setdiff(names(df), c(input$outcome, input$treatment))

    aipw_result <- AIPW$new(
      Y = df[[input$outcome]],
      A = df[[input$treatment]],
      X = df[covariates],
      Q.SL.library = c("SL.glm"),
      g.SL.library = c("SL.glm"),
      k_split = 2
    )
    aipw_result$fit()
    aipw_result$summary()$est
  })

  output$ate_output <- renderPrint({
    req(estimate_ate())
    cat("Estimated ATE:", estimate_ate())
  })
}

shinyApp(ui, server)
