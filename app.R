library(shiny)
library(ggplot2)

# Define UI
ui <- fluidPage(
  titlePanel("Advanced BED Calculator"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("m", "Number of dose fractions (m):", 30, min = 1),
      numericInput("d", "Dose/fraction (d) in Gy:", 2, min = 0),
      numericInput("T", "Total treatment time (T) in days:", 42, min = 1),
      hr(),
      numericInput("alpha", "Linear tumor cell killing parameter (α) in 1/Gy:", 0.2, min = 0, step = 0.01),
      numericInput("r", "Alpha/beta ratio (α/β) in Gy:", 10, min = 0, step = 0.1),
      numericInput("L_DI", "Accelerated repopulation rate (DI model):", 0.2, min = 0, step = 0.01),
      numericInput("L_DD", "Max repopulation rate (DD model):", 0.5, min = 0, step = 0.01),
      numericInput("Tk", "Time to accelerated repopulation (days):", 28, min = 0, step = 1)
    ),
    
    mainPanel(
      h3("Results:"),
      verbatimTextOutput("bedResults")
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  # Reactive expression for calculations
  bed_calcs <- reactive({
    m <- input$m
    d <- input$d
    T <- input$T
    alpha <- input$alpha
    r <- input$r
    L_DI <- input$L_DI
    L_DD <- input$L_DD
    Tk <- input$Tk
    
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
