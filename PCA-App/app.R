# ---------------------------------------
# Applet for principal component analysis
# ---------------------------------------
# 1. Simulate a 2D sample from a multivariate normal
# 2. Center it
# 3. Ask for an angle theta
# 4. Draw projected points on the line directed by that angle
#    and on the orthogonal
# 5. Print variances of projected points (inertia)
# 6. Maximize manually, and compare with theory: 
#    eigenvalue of the empirical covariance matrix
#    Observe that the two variances sum to the total variance.
# -----------------------------------------
# Author: O. Roustant, INSA Toulouse (2020)
#         roustant@insa-toulouse.fr
# -----------------------------------------

library(shiny)

ui <- fluidPage(
  titlePanel("Illustration of Principal Component Analysis (PCA)"),
  sidebarLayout(
    sidebarPanel(
      sliderInput(inputId = "n",
              label = "Sample size",
              value = 20, min = 10, max = 100, step = 10),
      sliderInput(inputId = "cor",
              label = "Correlation value",
              value = 0.8, min = - 1, max = 1, step = 0.1),
      sliderInput(inputId = "varRatio",
              label = "Ratio of variance (first input vs second one)",
              value = 2, min = 0.2, max = 5),
      hr(),
      sliderInput(inputId = "angle",
              label = "Can you find the angle for which the variance of projections is maximal?",
              value = 0, min = -90, max = 90),
      hr(),
      radioButtons(inputId = "solFlag",
               label = "Solution: eigen decomposition 
               of the empirical covariance matrix",
               choiceNames = list("Let me find the optimal angle!",
                                  "Show solution"),
               choiceValues = list(FALSE, TRUE))
    ),
    mainPanel(
      h4("The applet visualizes in 2D the solution of the PCA problem, 
          given by the diagonalization of the empirical covariance matrix."),
      h4("A sample is drawn from a bivariate Normal distribution.
          You have to find by yourself, only by a visual guess, 
          the first principal component (not a big deal here!).
          Then (2D case), the second one is given by the orthogonal axis."),
      h4("The two plots represent the two first principal axis, 
         and the projected data. 
         The numbers indicated are the proportion of variance of the projections,
         relatively to the total variance 
         (= sum of the variances of projected data onto the two axis)."),
      hr(),
      plotOutput(outputId = "PCA2Dplot"),
      verbatimTextOutput(outputId = "stats")
    )
  )
)

# Define server logic required to plot various variables against mpg
server <- function(input, output){
  source("PCA_util.R")
  library(MASS)
  data <- reactive(init(n = input$n, 
                        var1 = input$varRatio,
                        var2 = 1,
                        cor12 = input$cor))
  output$PCA2Dplot <- renderPlot(
    PCAturningLine(data(), 
                   solFlag = input$solFlag,
                   angle = input$angle)
  )
  output$stats <- renderPrint(
    if (input$solFlag) {
      SigmaHat <- cov(data()) * (input$n - 1) / input$n 
      E <- eigen(SigmaHat, symmetric = TRUE)
      cat("\nEigenvalues:\n")
      print(E$values)
      cat("\nCorresponding percentage of their sum:\n")
      print(round(100 * E$values / sum(E$values), 2))
      cat("\nEigenvectors:\n")
      print(E$vectors)
      angle <- atan(E$vectors[2 , 1]/E$vectors[1, 1]) * 180 / pi
      cat("\nCorresponding angle (in degrees): ", round(angle, 2))
    }
  )
}


shinyApp(ui = ui, server = server)
