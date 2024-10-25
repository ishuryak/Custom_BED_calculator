

library(shiny)
library(ggplot2)

# Define UI
ui <- fluidPage(
  titlePanel("BED Calculator"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("m", "Number of dose fractions (m):", 30, min = 1),
      numericInput("d", "Dose/fraction (d) in Gy:", 2, min = 0),
      numericInput("T", "Total treatment time (T) in days:", 42, min = 1)
    ),
    
    mainPanel(
      h3("Results:"),
      verbatimTextOutput("bedResults")
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  # Parameters
  alpha <- 0.2
  r <- 10
  L_DI <- 0.2
  L_DD <- 0.5
  Tk <- 28
  
  # Reactive expression for calculations
  bed_calcs <- reactive({
    m <- input$m
    d <- input$d
    T <- input$T
    
    BED_DI <- m*d*(1 + d/r) - L_DI*log(1 + exp(3*T - 3*Tk))/(3*alpha)
    
    BED_DD <- (r*L_DD*(exp(-m*alpha*d*(d + r)/(T*r)) - 1)*log(1 + exp(3*(47*m*d^2 + 47*d*m*r 
                                                                         + 28*r*log((1 + exp(141 - 3*Tk))^(-L_DI/((exp(-(84*alpha)/47) - 1)*L_DD)) - 1) 
                                                                         - 3948*r)*T/(47*m*d*(d + r)))) + 3*m*alpha*d*(d + r))/(3*r*alpha)
    
    list(BED_DI = BED_DI, BED_DD = BED_DD)
  })
  
  # Output
  output$bedResults <- renderPrint({
    results <- bed_calcs()
    cat("BED_DI:", round(results$BED_DI, 2), "\n")
    cat("BED_DD:", round(results$BED_DD, 2))
  })
}

shinyApp(ui = ui, server = server)
