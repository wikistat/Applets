# ---------------------------------------
# Applet for linear discriminant analysis
# ---------------------------------------
# -----------------------------------------
# Author: O. Roustant, INSA Toulouse (2020)
#         roustant@insa-toulouse.fr
# -----------------------------------------

library(shiny)

ui <- fluidPage(
  titlePanel("Illustration of Linear Discriminant Analysis (LDA)"),
  sidebarLayout(
    sidebarPanel(
      sliderInput(inputId = "n",
              label = "Class size",
              value = 20, min = 10, max = 100, step = 10),
      sliderInput(inputId = "nClass",
                  label = "Number of classes",
                  value = 3, min = 2, max = 8, step = 1),
      sliderInput(inputId = "rhoW",
              label = "Within-class correlation value",
              value = 0.8, min = - 1, max = 1, step = 0.1),
      sliderInput(inputId = "varW",
              label = "Within-class variance",
              value = 0.2, min = 0.1, max = 1, step = 0.1),
      hr(),
      radioButtons(inputId = "pred",
               label = "Exploration / Prediction",
               choiceNames = list("Exploration",
                                  "Show prediction frontiers"),
               choiceValues = list(FALSE, TRUE))
    ),
    mainPanel(
      h4("The applet shows the result of LDA, in a exploration / prediction perspective."),
      h4("Samples of equal size are drawn from a bivariate Normal distribution 
      with the same covariance matrix (within-class), 
      and different means drawn from a bivariate Normal distribution
      with identity covariance matrix (between-class)."),
      h4("For exploration, we draw the principal components 
         of the centroids with Mahalanobis metric,
         both for original data (left) and sphered data (right), 
         obtained by reducing with the estimated within-class covariance matrix."),
      h4("Prediction frontiers can be shown. When classes have the same size, 
         the optimal rule (Bayes classifier) is to compute, for sphered data, 
         the hyperplane bisectors of the sphered centroids (in 2D: Voronoi tesselation)."),
      hr(),
      plotOutput(outputId = "LDAplot")
    )
  )
)

# Define server logic required to plot various variables against mpg
server <- function(input, output){
  source("LDA_util.R")
  library(MASS)
  library(deldir)
  output$LDAplot <- renderPlot(
    LDAplot(nClass = input$nClass, 
            n = input$n, 
            rhoW = input$rhoW, 
            varW = input$varW,
            pred = input$pred)
  )
}


shinyApp(ui = ui, server = server)
