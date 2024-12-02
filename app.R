library(shiny)
library(kableExtra)
library(gridExtra)
library(ggplot2)
library(shinythemes)
library(shinyWidgets)
library(rsconnect)
source("required_packages.R")
source("Functions/support_functions.R")
source("Functions/initialization_functions.R")
source("Functions/E-step.R")
source("Functions/M-step.R")
source("Functions/log_likelihood.R")
source("Functions/EM-algorithm.R")
source("Functions/BIC_calculation.R")
source("Functions/table_BIC_comparison.R")

pedix <- c("\u2081","\u2082","\u2083","\u2084","\u2085","\u2086","\u2087","\u2088","\u2089")

b_rgb <- col2rgb("blue")
dp_rgb <- col2rgb("deeppink")
or_rgb <- col2rgb("orange")
g_rgb <- col2rgb("green")
r_rgb <- col2rgb("red")
p_rgb <- col2rgb("purple")

my_colors <- c("rgba(0, 0, 255, 0.6)", "rgba(0, 255, 0, 0.6)",
               "rgba(255, 20, 147, 0.6)", "rgba(255, 165, 0, 0.6)", "rgba(255, 0, 0, 0.6)",
               "rgba(160, 32, 240, 0.6)")

#my_colors <- c(alpha("blue",0.4),my_colors[2],my_colors[3],
               #alpha("orange",0.4))
border_string <- "border-bottom: 2px solid;border-top: 2px solid;border-left: 2px solid;border-right: 2px solid;"

# Define UI
ui <- fluidPage(theme = shinytheme("cerulean"),
                setBackgroundColor("white"), 
  navbarPage(div("ML-CWMd Model", style = "width: 240px; font-size: 28px; height: 22px; font-family: cursive;"),
             tabPanel(div("Upload Data  & Run Algorithm", style = "width: 265px; font-size: 20px; height: 22px;"), 
                      fluidPage(
                        sidebarLayout(
                          sidebarPanel(
                            fileInput("file", div(class = "content2", tags$ul(
                              tags$li("Choose a CSV file"),
                            )), accept = ".csv"),
                            tags$style("
                             .btn-file {  
                     background-color:#42B3FF; 
                       border-color: #42B3FF; 
                      }

             .progress-bar {
             background-color: #42B3FF;
             }

             "),
                            tags$style(type="text/css",
                                       ".shiny-output-error { visibility: hidden; }",
                                       ".shiny-output-error:before { visibility: hidden; }"
                            ),
                            # p(HTML("<b>After uploading the CSV file, you have the choice among these options:</b>")),
                            # tags$ul(
                            #   style = "padding-left: 10px;",
                            #   tags$li(
                            #     HTML("<span style='background-color: #3498DB; color: white; padding: 2px; border-radius: 3px;font-family: cursive; border: 0px solid black; background: linear-gradient(156.8deg, rgb(30, 144, 231) 27.1%, rgb(67, 101, 225) 77.8%);'></b>Continuous</b></span> : select which variables to model as continuous variables")
                            #   ),
                            #   tags$br(),
                            #   tags$li(
                            # HTML("<span style='background-color: #3498DB; color: white; padding: 2px; border-radius: 3px;font-family: cursive; background: linear-gradient(156.8deg, rgb(30, 144, 231) 27.1%, rgb(67, 101, 225) 77.8%);'></b>Categorical</b></span> : select which variables to model as categorical variables")),
                            #   tags$br(),
                            #   tags$li(HTML("<span style='background-color: #3498DB; color: white; padding: 2px; border-radius: 3px;font-family: cursive;background: linear-gradient(156.8deg, rgb(30, 144, 231) 27.1%, rgb(67, 101, 225) 77.8%);'></b>Binary Dependent</b></span> : select which variables to model as binary dependent variables")),
                            #   tags$br(),
                            #   tags$li(HTML("<span style='background-color: #3498DB; color: white; padding: 2px; border-radius: 3px;font-family: cursive;background: linear-gradient(156.8deg, rgb(30, 144, 231) 27.1%, rgb(67, 101, 225) 77.8%);'></b>Fixed Effects</b></span> : select which variables to include as fixed-effects in the regression part of the model")),
                            #   tags$br(),
                            #   tags$li(HTML("<span style='background-color: #3498DB; color: white; padding: 2px; border-radius: 3px;font-family: cursive;background: linear-gradient(156.8deg, rgb(30, 144, 231) 27.1%, rgb(67, 101, 225) 77.8%);'></b>Response Variable</b></span> : select the binary response variable of the regression part")),
                            #   tags$br(),
                            #   tags$li(HTML("<span style='background-color: #3498DB; color: white; padding: 2px; border-radius: 3px;font-family: cursive;background: linear-gradient(156.8deg, rgb(30, 144, 231) 27.1%, rgb(67, 101, 225) 77.8%);'></b>Grouping Level</b></span> : select the variable containing the group levels"))
                            # ) ,
                            fluidRow(
                              column(12, uiOutput("clusters")),  
                            ),
                            fluidRow(
                              column(12, uiOutput("explanation")),  
                            ),
                            fluidRow(
                              column(12, uiOutput("init")),  
                            ),
                            fluidRow(
                              column(12, uiOutput("nchoice")),  
                            ),
                      fluidRow(
                              column(3, uiOutput("columnSelection1")),  # First column selection
                              column(3, uiOutput("columnSelection2")),
                              column(3, uiOutput("columnSelection3")),
                              column(3, uiOutput("columnSelection4")),
                            ),
                            fluidRow(  
                              column(6, uiOutput("singleColumnSelection")),
                              column(6, uiOutput("singleColumnSelection1"))# Fourth column selection
                            ),
                            actionButton("runAlgorithm", "Run Algorithm", icon("paper-plane"), 
                                         style="color:#42B3FF ; background-color:#42B3FF ; border-color: #42B3FF "),
                            textOutput("uploadStatus"),
                      textOutput("catOutput"), # Output for cat() result
                            
                          width = 6),
                          mainPanel(
                            tags$head(
                              tags$link(rel="stylesheet", 
                                        href="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/katex.min.css", 
                                        integrity="sha384-dbVIfZGuN1Yq7/1Ocstc1lUEm+AT+/rCkibIcC/OmWo5f0EA48Vf8CytHzGrSwbQ",
                                        crossorigin="anonymous"),
                              HTML('<script defer src="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/katex.min.js" integrity="sha384-2BKqo+exmr9su6dir+qCw08N2ZKRucY4PrGQPPWU1A7FtlCGjmEGFqXCv5nyM5Ij" crossorigin="anonymous"></script>'),
                              HTML('<script defer src="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/contrib/auto-render.min.js" integrity="sha384-kWPLUVMOks5AQFrykwIup5lo0m3iMkkHrD0uJ4H5cjeGihAutqP0yW0J6dpFiVkI" crossorigin="anonymous"></script>'),
                              HTML('
    <script>
      document.addEventListener("DOMContentLoaded", function(){
        renderMathInElement(document.body, {
          delimiters: [{left: "$", right: "$", display: false}]
        });
      })
    </script>'),
          tags$style(HTML('
          .subtitl {
        color: black;
        margin-top: 20px;
        font-size: 20px;
        line-height: 1.6;
        text-align: center;
      }
          .titl {
        color: black;
        margin-top: 20px;
        font-size: 26px;
        line-height: 1.6;
        text-align: center;
          }
      .titl2 {
        color: black;
        margin-top: 20px;
        font-size: 26px;
        margin-left: 16px;
        line-height: 1.6;
      }
        .equation {
        color: black;
        margin-top: 20px;
        font-size: 16px;
        line-height: 1.6;
        text-align: center;
      }
      .content {
        color: black;
        margin-top: 20px;
        font-size: 16px;
        line-height: 1.6;
      }
      .content ul {
        list-style: none;
        margin-left: 0;
        padding-left: 20px;
      }
      .content ul li {
        margin-bottom: 10px;
        position: relative;
        padding-left: 30px; /* Increased padding for more space */
        color: black; /* Ensures the text is black */
      }
      .content ul li::before {
        content: "➤";
        color: #42B3FF;
        font-size: 18px;
        position: absolute;
        left: 0;
        top: 0;
      }
      
      .content2 {
        color: black;
        margin-top: 3px;
        font-size: 14px;
        line-height: 0.6;
        font-family: arial;
      }
      .content2 ul {
        list-style: none;
        margin-left: 0;
        padding-left: 0px;
      }
      .content2 ul li {
        margin-bottom: 5px;
        position: relative;
        padding-left: 19px; 
        color: #666666; 
      }
      .content2 ul li::before {
        content: "➤";
        color: #42B3FF;
        font-size: 14px;
        position: absolute;
        left: 0;
        top: 0;
      }
      .note {
      font-size: 12px;
      font-family: cursive;
        background-color: #42B3FF; 
        color: white; 
        padding-top: 10px; 
        border-radius: 5px; 
        margin-bottom: 0px; 
        text-align: center; 
        border: 3px solid #42B3FF; 
        display: flex;
        align-items: center;
        justify-content: center;
        width: 612px; 
         height: 19px; 
        margin-left: 22px;
      }
      .note2 {
      font-size: 12px;
      font-family: cursive;
        background-color: #42B3FF; 
        color: white; 
        padding-top: 10px; 
        border-radius: 5px; 
        margin-bottom: 0px; 
        border: 3px solid #42B3FF; 
        display: flex;
        align-items: center;
        width: 434px; 
         height: 19px; 
        margin-left: 22px;
      }
      .note3 {
      font-size: 12px;
      font-family: cursive;
        background-color: #42B3FF; 
        color: white; 
        padding-top: 10px; 
        border-radius: 5px; 
        margin-bottom: 0px; 
        border: 3px solid #42B3FF; 
        display: flex;
        align-items: center;
        width: 308px; 
         height: 19px; 
        margin-left: 22px;
      }
    ') )
                            ),
                            titlePanel(
                              div(class = "titl",HTML("<span style='font-family: cursive; color: white; background-color: #42B3FF; padding: 3px; border-radius: 5px; box-shadow: 0px 4px 6px rgba(0, 0, 0, 0.1);'>Model Description"))),
          div(class = "content",
              p("R implementation of ML-CWMd, a model rooted in the cluster-weighted models domain. It's designed for scenarios featuring:", style = "color: black;"),
              tags$ul(
                tags$li("Hierarchical data"),
                tags$li("Binary response"),
                tags$li("Assumed presence of latent clusters among observations"),
                tags$li("Multiple types of covariates: continuous, independent categorical, and dichotomous dependent covariates are all supported by the model.")
              )),
          div(class = "subtitl",
              p(HTML("<span style='font-family: cursive; color: white; background-color: #42B3FF; padding: 3px; border-radius: 5px; box-shadow: 0px 4px 6px rgba(0, 0, 0, 0.1);font-size: 22px;'>EM-Algorithm for Model Estimation"))),
          div(class = "titl",
                            img(src = "img.png", height="100%", width="100%", align = "center")),
                            p(""),
          div(class = "subtitl",
                            p(HTML("<span style='font-family: cursive; color: white; background-color: #42B3FF; padding: 3px; border-radius: 5px; box-shadow: 0px 4px 6px rgba(0, 0, 0, 0.1);font-size: 22px;'>Variables"))),
                            helpText(
                              div(class = "content",
                              tags$ul(
                                tags$li("$U$ represents $p$-dimensional vector of continuous variables"),
                                tags$li("$V$ represents $q$-dimensional vector of categorical covariates"),
                                tags$li("$D$ represents $h$-dimensional vector of dichotomous covariates which may possess some degree of dependence"),
                                tags$li("$X = (U,V,D)$ and $Y$ (response variable) ∈ $\\mathbb{R}^{(p + q + h)}$ × {0,1} defined in a finite space $\\boldsymbol{\\Omega}$ which is assumed to be divided into $C$ clusters denoted as $\\boldsymbol{\\Omega}_1,\\dots, \\boldsymbol{\\Omega}_C$"),
                                tags$li("Two-levels hierarchy, first-level observations $i$, for $i = 1, \\dots, n_j$, are nested within groups $j$, for $j = 1, \\dots, J$")
                              ), style = "color: black;"
                            )),
                            helpText(
                              div(class = "subtitl",   p(p(HTML("<span style='font-family: cursive; color: white; background-color: #42B3FF; padding: 3px; border-radius: 5px; box-shadow: 0px 4px 6px rgba(0, 0, 0, 0.1);font-size: 22px;'>Joint probability across the clusters")))),
                              div(class = "equation", p("$ p((\\boldsymbol{x},\\boldsymbol{y})| \\boldsymbol{\\theta}) = \\sum_{c= 1}^{C}
p(\\mathbf{y}|\\mathbf{x},\\boldsymbol{\\xi}_{c}) \\phi( \\mathbf{u}|
                       \\boldsymbol{\\mu}_{c}, \\boldsymbol{\\Sigma}_c ) \\psi(\\mathbf{v}|\\boldsymbol{\\lambda}_{c}) \\zeta(\\boldsymbol{d}| \\boldsymbol{\\Gamma}_c,\\boldsymbol{\\nu}_c)  w_{c} $"),
                             style = "color: black;")),
          div(class = "content",
                            helpText(tags$ul(
                              tags$style(HTML(".shiny-text-output { color: black; }")),
                              tags$li("$\\boldsymbol{\\theta}$ → vector containing all the parameters of the model"),
                              tags$li("$p(\\mathbf{y}|\\mathbf{x},\\bm{\\xi}_{c})$ → represents a multilevel logistic regression model where $\\boldsymbol{\\xi}_c$ are the cluster-specific parameters for both fixed and random effects"),
                              tags$li("$w_c$ → indicates the proportion of observations within cluster c"),
                              tags$li("$\\phi(\\cdot|\\boldsymbol{\\mu}_c,\\boldsymbol{\\Sigma}_c)$ → multivariate normal density with cluster-wise different mean vectors $\\boldsymbol{\\mu}_c$ and covariance matrices $\\boldsymbol{\\Sigma}_c$"),
                              tags$li("$\\psi(\\cdot|\\boldsymbol{\\lambda}_c)$ → independent multinomial distributions with cluster-wise different parameter vectors $\\boldsymbol{\\lambda}_c$"),
                              tags$li("$\\zeta(\\cdot| \\boldsymbol{\\Gamma}_c,\\boldsymbol{\\nu}_c)$ → Ising model with cluster-wise different threshold vectors $\\boldsymbol{\\nu}_c$ and interaction matrices $\\boldsymbol{\\Gamma}_c$")
                            ), style = "color: black;")),
                            
                          width = 6)
                        )
                      )
             ),
             tabPanel(div("Model Selection", style = "width: 145px; font-size: 20px; height: 22px;"),
                      column(12,uiOutput("biccomp"))),
             tabPanel(div("Results Visualization", style = "width: 185px; font-size: 20px; height: 22px;"),
                      fluidPage(
                        sidebarLayout(
                          sidebarPanel(
                            h4("Model Parameters"),
                            selectInput("choice", "Select Output:", choices = c(
                              "Continuous Covariates Parameters",
                              "Categorical Covariates Parameters",
                              "Random Effects",
                              "Fixed Effects",
                              "Ising Model"
                            )), width = 2
                          ),
                          mainPanel(
                            textOutput("resultText"),
                            tableOutput("outputTable"),
                            conditionalPanel(
                              condition = "input.choice == 'Ising Model'",
                              column(12,uiOutput("NuuTables")),
                              column(12,uiOutput("IntPlot")),
                              # fluidRow(
                              #   column(4, plotOutput("plot1")), 
                              #   column(4, plotOutput("plot2")),  
                              #   column(4, plotOutput("plot3"))   
                              # )
                            ),
                            conditionalPanel(
                              condition = "input.choice == 'Random Effects'",
                              column(12,uiOutput("RanefTables")),
                              column(12,uiOutput("RanefPlot")),
                              # fluidRow(
                              #   column(4, plotOutput("plotre1")),  
                              #   column(4, plotOutput("plotre2")), 
                              #   column(4, plotOutput("plotre3"))   
                              # )
                            ),
                            conditionalPanel(
                              condition = "input.choice == 'Continuous Covariates Parameters'",
                              uiOutput("ContTables")
                            ),
                            conditionalPanel(
                              condition = "input.choice == 'Fixed Effects'",
                              uiOutput("fixedEffectsTables")
                            ),
                            conditionalPanel(
                              condition = "input.choice == 'Categorical Covariates Parameters'",
                              uiOutput("CatTables")
                            ), width = 10
                            
                            # Add other conditional panels for different model parameters as before
                          )
                        )
                      )
                      
             ),
             tabPanel(div("Model Predictions", style = "width: 160px; font-size: 20px; height: 22px;"),
             fluidPage(
               sidebarLayout(
                 sidebarPanel(
                   fileInput("file_pred", HTML("<span style='font-family: cursive;'>Choose a CSV file containing new observations"), accept = ".csv"),
                   tags$style("
                             .btn-file {
                     background-color:#42B3FF;
                       border-color: #42B3FF;
                      }

             .progress-bar {
             background-color: #42B3FF;
             }

             "),
                   tags$style(type="text/css",
                              ".shiny-output-error { visibility: hidden; }",
                              ".shiny-output-error:before { visibility: hidden; }"
                   ),
                   actionButton("runPrediction", "Run Prediction", icon("paper-plane"),
                                style="color:#3498DB ; background-color:#3498DB ; border-color: #3498DB "),
                   textOutput("uploadStatusPred"),
                   textOutput("catOutputPred"), # Output for cat() result

                   width = 4),
                 mainPanel(
                   titlePanel(div(class = "titl2",HTML("<span style='font-family: cursive; color: white; background-color: #42B3FF; padding: 3px; border-radius: 5px; box-shadow: 0px 4px 6px rgba(0, 0, 0, 0.1);'>Model Predictions"))),
                   column(8, uiOutput("tablepred")),
                   width = 8)
               )
             ) )
  )
)



server <- function(input, output, session) {
 
  
  data <- reactive({
    req(input$file)
    read.csv(input$file$datapath)
  })
  
  # Display upload status
  output$uploadStatus <- renderText({
    if (is.null(input$file))
      return("Please upload a CSV file.")
    else
      return("File uploaded successfully.")
  })
  
  data_pred <- reactive({
    req(input$file_pred)
    read.csv(input$file_pred$datapath)
  })
  
  output$uploadStatusPred <- renderText({
    if (is.null(input$file_pred))
      return("Please upload a CSV file.")
    else
      return("File uploaded successfully.")
  })
  
  observeEvent(input$file, {
    output$clusters <- renderUI({
      if (is.null(input$file)) return(NULL)
      checkboxGroupInput("ClusterChoice",  div(class = "content2", tags$ul(
        tags$li("Should the number of clusters (C) be fixed or vary on a range?"),
      )) , choices = c("Fixed", "Range"), selected = NULL)
    })
  })
  
  observe({
    req(input$file)
    
    if (!is.null(input$file)) {
      clusters <- req(input$ClusterChoice)
      
      if ("Fixed" %in% clusters) {
        output$explanation <- renderUI({
          textAreaInput("userInput", div(class = "content2", tags$ul(
            tags$li("Enter the number of clusters (e.g., 3):"),
          )), value = "", width = "700px", height = "40px")
        })
      } else if ("Range" %in% clusters) {
        output$explanation <- renderUI({
          textAreaInput("userInput", helpText(div(class = "content2", tags$ul(
            tags$li("Specify the range (e.g., 2-4) to run the model across a series of values (e.g., C = 2, 3, 4)"),
          )), div(class = "note", 
            p("In the 'Model Selection' section, you'll find the best results based on BIC within the specified range"),
          ) ) , value = "",
                        width = "700px", height = "40px")
        })
      }
    }
  })
  
  observeEvent(input$file, {
    output$init <- renderUI({
      if (is.null(input$file)) return(NULL)
      checkboxGroupInput("InitChoice", div(class = "content2", tags$ul(
        tags$li("Would you like to initialize the algorithm with random initialization or k-means initialization?"),
      )) , 
                         choices = c("Random", "k-means"), selected = NULL)
    })
  })
  
  observeEvent(input$file, {
    output$nchoice <- renderUI({
      if (is.null(input$file)) return(NULL)
      textAreaInput("userN", helpText(div(class = "content2", tags$ul(
        tags$li("Specify the number of initializations M to perform (e.g., 20)"),
      )), div(class = "note2", 
              p("For random initialization, 
                                   it's recommended to set M between 20 and 50"),
      ),
      p(""),
      div(class = "note3", 
          p("For k-means initialization, set M between 5 and 10"),
      )), value = "",
                    width = "700px", height = "40px")
    })
  })
  
  user_text <- reactive({
    input$userInput
  })
  
  user_init <- reactive({
    input$InitChoice
  })
  
  
  user_n <- reactive({
    input$userN
  })
  
  
  # Dynamic UI for selecting columns after data upload - Column Selection 1
  observeEvent(input$file, {
    output$columnSelection1 <- renderUI({
      if (is.null(input$file)) return(NULL)
      checkboxGroupInput("selectedColumns1", div(class = "content2", tags$ul(
        tags$li("Continuous"),
      )), 
                         choices = colnames(data()), selected = NULL)
    })
  })
  
  # Dynamic UI for selecting columns after data upload - Column Selection 2
  observeEvent(input$file, {
    output$columnSelection2 <- renderUI({
      if (is.null(input$file)) return(NULL)
      checkboxGroupInput("selectedColumns2",  div(class = "content2", tags$ul(
        tags$li("Categorical"),
      )),
                         choices = colnames(data()), selected = NULL)
    })
  })
  
  # Dynamic UI for selecting columns after data upload - Column Selection 3
  observeEvent(input$file, {
    output$columnSelection3 <- renderUI({
      if (is.null(input$file)) return(NULL)
      checkboxGroupInput("selectedColumns3",  div(class = "content2", tags$ul(
        tags$li("Binary Dependent"),
      )), 
                         choices = colnames(data()), selected = NULL)
    })
  })
  
  # Dynamic UI for selecting columns after data upload - Column Selection 4
  observeEvent(input$file, {
    output$columnSelection4 <- renderUI({
      if (is.null(input$file)) return(NULL)
      checkboxGroupInput("selectedColumns4", div(class = "content2", tags$ul(
        tags$li("Fixed Effects"),
      )), 
                         choices = colnames(data()), selected = NULL)
    })
  })
  
  # Dynamic UI for selecting a single column after data upload
  observeEvent(input$file, {
    output$singleColumnSelection <- renderUI({
      if (is.null(input$file)) return(NULL)
      selectInput("singleColumn",  div(class = "content2", tags$ul(
        tags$li("Response Variable"),
      )), 
                  choices = colnames(data()), selected = colnames(data())[1])
    })
  })
  
  # Dynamic UI for selecting a single column after data upload
  observeEvent(input$file, {
    output$singleColumnSelection1 <- renderUI({
      if (is.null(input$file)) return(NULL)
      selectInput("singleColumn1", div(class = "content2", tags$ul(
        tags$li("Grouping Level"),
      )), 
                  choices = colnames(data()), selected = colnames(data())[1])
    })
  })
  
  
  selected_data1 <- reactive({
    if (!is.null(input$file)) {
      data() %>% select(input$selectedColumns1)
    }
  })
  
  #--------------------------------RUN ALGORITHM--------------------------------
  cat_output <- reactiveVal(NULL)
  
  fitting <- NULL
  U <-  NULL
  V <-  NULL
  V1 <-  NULL
  V2 <-  NULL
  D <-  NULL
  Y <-  NULL
  formula <-  NULL
  C <-  NULL
  mydata <-  NULL
  
  log <- reactiveVal("") 
  
  observeEvent(input$runAlgorithm, {
    
    Calculate_BIC <- function(formula, data, U = NULL, V = NULL, D = NULL, log_l, C, fixed_intercept = FALSE){
      N = dim(data)[1]
      k_weights = C - 1
      
      if(!is.null(U)){
        U <- data.frame(U)
        p = dim(U)[2]
        k_cont = p*(p+3)/2
      } else{
        k_cont <- 0
      }
      
      if(!is.null(V)){
        k_r = unlist(purrr::map(V, ~dim(.)[2]))
        k_cat = sum(k_r)
      } else{
        k_cat <- 0
      }
      
      if(!is.null(D)){
        l = dim(D)[2]
        k_dich = l*(l+1)/2
      } else{
        k_dich <- 0
      }
      
      k_reg = length(attr(terms(formula), "term.labels")) + ifelse(fixed_intercept == TRUE, 1, 0)
      
      k = C*(k_cont + k_cat + k_dich + k_reg) + k_weights
      
      BIC <- -2*log_l + log(N) * k
      return(BIC)
      
    }
    
    fitMLCWMd <- function(data, formula, C, cont_var = NULL, cat_var = NULL, dich_dep_var = NULL, 
                          init, cluster = NULL, mstep_isingfit = FALSE, 
                          fixed_intercept = TRUE,
                          max_it = 30, tol = 1e-6){
      
      # Process variables
      if(length(cont_var) != 0){
        missing_cont_var <- setdiff(cont_var, colnames(data))
        if(length(missing_cont_var) == 0){
          U <- data[ , cont_var]
        } else{
          return("At least one continuous variable provided is not in the data")
        }
      } else{
        U <- NULL
      }
      
      if(length(cat_var) != 0){
        missing_cat_var <- setdiff(cat_var, colnames(data))
        if(length(missing_cat_var) == 0 ){
          cat_names_list <- list()
          for(i in 1:length(cat_var)){
            #cat_names_list[[i]] <- paste0(cat_var[i], 1:length( levels( as.factor(data[, cat_var[i]]) ) ) )
            cat_names_list[[i]] <- paste0(cat_var[i], ".", levels( as.factor(data[, cat_var[i]])  ) )
          }
          cat_data <- as.data.frame(data[ , cat_var])
          colnames(cat_data) <- cat_var
          cat_data <- as.data.frame(unclass(cat_data),stringsAsFactors = TRUE)
          onehot_cat_data <- cat_to_onehot(cat_data, cat_var, cat_names_list)
          V <- list()
          for(i in 1:length(cat_var)){
            V[[i]] <- onehot_cat_data[, grep(paste0("^", cat_var[i]), names(onehot_cat_data))]
          }
        } else{
          return("At least one categorical variable provided is not in the data")
        }
      } else{
        V <- NULL
      } 
      
      if(length(dich_dep_var) != 0){
        missing_dich_dep_var <- setdiff(dich_dep_var, colnames(data))
        if(length(missing_dich_dep_var) == 0){
          D <- data[ , dich_dep_var]
        } else{
          return("At least one dichotomous dependent variable provided is not in the data")
        }
      } else{
        D <- NULL
      }
      
      Y <- data[ , as.character(formula[[2]])]
      data <- as.data.frame(unclass(data),stringsAsFactors = TRUE)
      #return(list(U,V,D))
      
      # Decide Initialization type
      if(init == "random"){
        params = initialize_params_random(data,U,V,D,formula,C)
      } else if (init == "kmeans"){
        params = initialize_params_k_means(data,U,V,D,formula,C)
      } else if (init == "manual"){
        params = initialize_params_manual(data,cluster,U,V,D,formula,C)
      } else{
        stop("Execution stopped. No existing initialization provided")
      }
      
      z <- E_step(Y,U,V,D,C,params)
      log_l <- c(0,log_like(Y,U,V,D,C,z,params))
      tol = tol
      max_it = max_it
      n_iter = 0
      while( ( abs(log_l[length(log_l)] - log_l[length(log_l)-1])  >  tol ) & ( n_iter < max_it ) ){
        # E-step
        z <- E_step(Y,U,V,D,C,params)
        # M-step
        if(mstep_isingfit == FALSE){
          params <- M_step_EI(formula,U,V,D,data,C,z,params)
        }
        else{
          params <- M_step_IF(formula,U,V,D,data,C,z,params)
        }
        #Saving log-likelihood
        new_log_l <- log_like(Y,U,V,D,C,z,params)
        log_l <- c(log_l,new_log_l)
        n_iter <- n_iter + 1
      }
      model_bic = Calculate_BIC(formula, data, U, V, D, log_l[length(log_l)], C, fixed_intercept = fixed_intercept)
      
      if(length(cont_var) != 0){
        dimnames(params$mu) <- list(NULL, cont_var, NULL)
        dimnames(params$sigma) <- list(cont_var, cont_var, NULL)
      } 
      
      if(length(cat_var) != 0){
        for(i in 1:length(cat_var)){
          dimnames(params$lambda[[i]]) <- list(NULL, cat_names_list[[i]], NULL)
        }
      }
      if(length(dich_dep_var) != 0){
        dimnames(params$thres) <- list(NULL, dich_dep_var, NULL)
        dimnames(params$int) <- list(dich_dep_var, dich_dep_var, NULL)
      }
      
      return(list(params = params,log_l = log_l, z = z, bic = model_bic,
                  data_spez = list(U = U, V = V, D = D)))
    }
    
    selected_continuous_covariates <- input$selectedColumns1
    selected_categorical_covariates <- input$selectedColumns2
    selected_binary_dependent_covariates <- input$selectedColumns3
    selected_fixed_effects <- input$selectedColumns4
    selected_response_variable <- input$singleColumn
    selected_grouping_levels <- input$singleColumn1
    
    mydata <- data()
    original_formula <- y ~ (1|level)
    formula <- as.formula(paste(selected_response_variable,
                                "~ (1|", selected_grouping_levels, ")"))
    formula_update <- paste(
      "~ . +",
      paste(selected_fixed_effects, collapse = " + ")
    )
    formula <- update.formula(formula, formula_update)
    
    if(nchar(user_text()) == 1){
      C <- as.numeric(user_text())
    }
    else{
      range_parts <- unlist(strsplit(user_text(), "-"))
      C_vec <- range_parts[1]:range_parts[2]
    }
    
    flag_ini <- if_else(user_init() == "Random", "random", "kmeans")
    flag_ini <<- flag_ini
    
    if(nchar(user_text()) == 1){
      N <- user_n()
      log_ll <- vector(mode="numeric", length = N)
      fitting_rec <- list()
      for(i in 1:N){
        tryCatch({
          fit <- fitMLCWMd(data = mydata, formula = formula,
                           C = C, 
                           cont_var = selected_continuous_covariates, 
                           cat_var = selected_categorical_covariates, 
                           dich_dep_var = selected_binary_dependent_covariates, 
                           init = flag_ini, 
                           mstep_isingfit = TRUE, 
                           fixed_intercept = TRUE,
                           max_it = 30, 
                           tol = 1e-6)
          log_ll[i] <- fit$log_l[length(fit$log_l)]
          fitting_rec[[i]] <- fit
        }, error = function(err) {
          cat(paste("Error occurred at index", i, "Error message:", conditionMessage(err), 
                    "\n", "Next Initialization...", "\n"))
        }
        )
      }
      
      # Select the result with the highest log-likelihood
      log_ll <- if_else(log_ll == 0, -Inf, log_ll)
      fit <- fitting_rec[[which.max(log_ll)]]
     
      final_log_l <- fit$log_l[length(fit$log_l)]
      Bbic <- fit$bic
      fitting <<- fit
      BIC_vec <- Bbic
    }
    else{
      saved_fitting <- list()
      BIC_vec <- vector(mode = "numeric", length = length(C_vec) )
      for(k in 1:length(C_vec)){
        fit_list <- list()
        N <- user_n()
        cat(paste("Starting running with C =", C_vec[k], "\n \n"))
        log_ll <- vector(mode="numeric", length = N)
        for(i in 1:N){
          tryCatch({
            fit <-  fitMLCWMd(data = mydata, formula = formula,
                              C = C_vec[k], 
                              cont_var = selected_continuous_covariates, 
                              cat_var = selected_categorical_covariates, 
                              dich_dep_var = selected_binary_dependent_covariates, 
                              init = flag_ini, 
                              mstep_isingfit = TRUE, 
                              fixed_intercept = TRUE,
                              max_it = 30, 
                              tol = 1e-6)
            log_ll[i] <- fit$log_l[length(fit$log_l)]
            fit_list[[i]] <- fit
          }, error = function(err) {
            cat(paste("Error occurred at index", i,  "for number of clusters C = ",  C_vec[k], "|","Error message:", conditionMessage(err), 
                      "\n", "Next Initialization...", "\n \n"))
          }
          )
        }
        log_ll <- if_else(log_ll == 0, -Inf, log_ll)
        fitting <- fit_list[[which.max(log_ll)]]
        
        BIC_vec[k] <- fitting$bic
        saved_fitting[[k]] <- fitting 
      }
      C <- C_vec[which.min(BIC_vec)]
      fitting <<- saved_fitting[[which.min(BIC_vec)]]
    }
    
    
    #fit <- list(params = params,log_l = log_l,z = z)
    #fitting <<- fit
    U <<- data.frame(fitting$data_spez$U)
    V <<- fitting$data_spez$V
    V1 <<- fitting$data_spez$V[[1]]
    D <<- fitting$data_spez$D
    C <<- C
    formula <<- formula
    mydata <<- mydata
    
    if(nchar(user_text()) > 1){
      BIC_vec <<- BIC_vec
      C_vec <<- C_vec
    }
    if(nchar(user_text()) == 1){
      BIC_vec <<- BIC_vec
      C_vec <<- C
    }
    
    selected_continuous_covariates <<- input$selectedColumns1
    selected_categorical_covariates <<- input$selectedColumns2
    selected_binary_dependent_covariates <<- input$selectedColumns3
    selected_fixed_effects <<- input$selectedColumns4
    selected_response_variable <<- input$singleColumn
    selected_grouping_levels <<- input$singleColumn1
    
    output_text <- capture.output({
      cat("Model Estimation COMPLETED. Go in Results Visualization Section to see it")
    })
    cat_output(paste(output_text, collapse = "\n"))
    
  })
  
  #----------------------------RUN PREDICTION-----------------------------------
  observeEvent(input$runPrediction, {
    
    # Model Prediction
    predict_MLCWMd <- function(fit, data, C, new_level = FALSE){
      if(!is.null(fit$params$mu)){
        U <- data[ , dimnames(fit$params$mu)[[2]] ]
        cont_term <- list()
        for(i in 1:C){
          cont_term[[i]] <- dmvnorm(U, mean = 
                                      fit$params$mu[,,i], sigma = fit$params$sigma[,,i],log = 1)
        }
      }  else{
        cont_term <- rep(0,C)
      }
      
      if(!is.null(fit$params$lambda)){
        V <- list()
        v <- length(fit$params$lambda)
        for(i in 1:v){
          temp_var <- sub("\\..*", "", dimnames(fit$params$lambda[[i]])[[2]])[1]
          temp_data <- as.data.frame(data[, temp_var])
          colnames(temp_data) <- temp_var
          V[[i]] <- cat_to_onehot(temp_data, temp_var, list(dimnames(fit$params$lambda[[i]])[[2]]))
        }
        multinom_term <- list()
        for(i in 1:C){
          multinom_term[[i]] <- all_mydmultinom_log(V,fit$params$lambda,v,i) 
        } 
      } else{
        multinom_term <- rep(0,C)
      }
      
      if(!is.null(fit$params$thres)){
        D <- data[ , dimnames(fit$params$thres)[[2]]]
        dich_term <- list()
        for(i in 1:C){
          dich_term[[i]] <- log(apply(D, 1, function(row) IsingStateProb(row, 
                                                                         fit$params$int[,,i], fit$params$thres[,,i], beta = 1, responses = c(0L, 1L))))
        }
      } else{
        dich_term <- rep(0,C)
      } 
      
      
      weights_ = matrix(0, nrow = dim(data)[1], ncol = C)
      for(i in 1:C){
        weights_[, i] = exp( log(fit$params$w[i]) + cont_term[[i]] +
                               + multinom_term[[i]] + dich_term[[i]] ) 
        
      }
      y_pred <- 0
      if(new_level == F){
        for(i in 1:C){
          y_pred <- y_pred + weights_[,i]*predict(fit$params$models[[i]], data, type = "response", allow.new.levels = F) 
          
        }
      }
      else{
        for(i in 1:C){
          y_pred <- y_pred +  weights_[,i]*predict(fit$params$models[[i]], data, type = "response", allow.new.levels = T) 
          
        }
      }
      
      y_pred <- y_pred / rowSums(weights_)
      
      return(y_pred)
    }
    
    
    mydata_pred <- data_pred()
    yy_pred <- predict_MLCWMd(fitting, mydata_pred, C, new_level = FALSE)
    yy_pred <<- yy_pred
    mydata_pred <<- mydata_pred
    
    if(!is.null(fitting$params$mu)){
      u <- length(names(fitting$params$mu[,,1]))
    } else{
      u <- 0
    }
    if(!is.null(fitting$params$thres)){
      d <- length(names(fitting$params$thres[,,1]))
    } else{
      d <- 0
    }
    if(!is.null(fitting$params$lambda)){
      v <- length(fitting$params$lambda)
    } else{
      v <- 0
    }
    output$tablepred <- renderText({
      N <- length(yy_pred)
      border_string <- "font-family: cursive; border-bottom: 2px solid;border-top: 2px solid;border-left: 2px solid;border-right: 2px solid;"
      ff = as.matrix( cbind("New" = 1:N, Prediction = round(yy_pred,2), mydata_pred ) )
      cn <- colnames(ff)
      colnames(ff) <- NULL
      ff %>% kbl(align = c("c", "c","c", "c","c", "c","c", "c","c", "c","c","c", "c","c", "c","c")) %>%
        kable_styling(font_size = 20)  %>%
        row_spec(seq(1,length(yy_pred)), extra_css = "border-bottom: 2px solid;") %>%
        column_spec(2, background = "#42B3FF") %>%
        column_spec(1, background = alpha("grey90",0.2)) %>%
        column_spec(1:(d + u + v + 3 + 1), extra_css = "border-left: 2px solid;", width = "100px") %>%
        column_spec(d + u + v + 3 + 1, extra_css = "border-right: 2px solid;") %>%
        add_header_above(cn, extra_css = border_string, color = "black", 
                         background = alpha("white",0.2))
    })
    
    
    
  })
  
  
  output$catOutput <- renderText({
    cat_output()
  })
  
  
  output$outputTable <- renderText({
    choice <- input$choice
    if (choice == "Random Effects") {
      output$RanefTables <- renderUI({
        ui_elements <- list()
        ui_elements <- append(ui_elements, list(
          HTML(paste0(
            "<div style='text-align: center; padding: 20px;'>",
            "<div style='border: 2px solid black; background-color: #FFAC1C; padding: 5px; display: inline-block; border-radius: 10px; box-shadow: 0px 0px #000000; background-image: linear-gradient(90deg, #b4d8ff 0%, #9dc2ff 35%, #c5e4ff 100%);'>",
            "<span style='font-size: 28px; font-family: cursive; color: #000000;'>Random Effects Table", "</span>",
            "</div>",
            "</div>"
          )),
          tableOutput("tableranef"),
          HTML("<br>"),
          HTML("<br>")
        ))
        ui_elements
      }) 
      
      output$RanefPlot <- renderUI({
        ui_elements <- list()
        ui_elements <- append(ui_elements, list(
          HTML(paste0(
            "<div style='text-align: center; padding: 20px;'>",
            "<div style='border: 2px solid black; background-color: #FFAC1C; padding: 5px; display: inline-block; border-radius: 10px; box-shadow: 0px 0px #000000; background-image: linear-gradient(90deg, #b4d8ff 0%, #9dc2ff 35%, #c5e4ff 100%);'>",
            "<span style='font-size: 28px; font-family: cursive; color: #000000;'>Random Effects Plots", "</span>",
            "</div>",
            "</div>"
          ))
        ))
        for(c in 1:C){
        ui_elements <- append(ui_elements, list(
          plotOutput(paste0("plotre", c)),
          HTML("<br>"),
          HTML("<br>")
          ))
        }
        
        ui_elements
      })  
      
    }
    
    
    else if (choice == "Continuous Covariates Parameters") {
      output$ContTables <- renderUI({
        if (is.null(selected_continuous_covariates)) return(NULL)
        
        uu = c("mu","sigma")
        uu1 = c("Means","Covariance Matrices")
        ui_elements <- list()
        for (u in 1:2) {
          ui_elements <- append(ui_elements, list(
            HTML(paste0(
              "<div style='text-align: center; padding: 20px;'>",
              "<div style='border: 2px solid black; background-color: #FFAC1C; padding: 5px; display: inline-block; border-radius: 10px; box-shadow: 0px 0px #000000; background-image: linear-gradient(90deg, #b4d8ff 0%, #9dc2ff 35%, #c5e4ff 100%);'>",
              "<span style='font-size: 28px; font-family: sans-serif; color: #000000;'>Continuous Variables ", uu1[u], "</span>",
              "</div>",
              "</div>"
            )),
            tableOutput(paste0("table", uu[u])),
            HTML("<br>"),
            HTML("<br>")
          ))
        }
        ui_elements
      })  
      
    }
    
    
    else if (choice == "Fixed Effects") {
    output$fixedEffectsTables <- renderUI({
      ui_elements <- list()
      for (c in 1:C) {
        ui_elements <- append(ui_elements, list(
          HTML(paste0(
            "<div style='text-align: center; padding: 20px;'>",
            "<div style='border: 2px solid black; background-color:", my_colors[c] ,"; padding: 5px; display: inline-block; border-radius: 10px; box-shadow: 0px 0px #000000;'>",
            "<span style='font-size: 28px; font-family: cursive; color: #000000;'>Fixed Effects Cluster ", c, "</span>",
            "</div>",
            "</div>"
          )),
          tableOutput(paste0("tablefix", c)),
          HTML("<br>"),
          HTML("<br>")
        ))
      }
      ui_elements
    })
    }
    
    else if (choice == "Categorical Covariates Parameters") {
      output$CatTables <- renderUI({
        if (is.null(selected_categorical_covariates)) return(NULL)
        
        ui_elements <- list()
        for (var in 1:length(V)) {
          ui_elements <- append(ui_elements, list(
            HTML(paste0(
              "<div style='text-align: center; padding: 20px;'>",
              "<div style='border: 2px solid black; background-color: #FFAC1C; padding: 5px; display: inline-block; border-radius: 10px; box-shadow: 0px 0px #000000; background-image: linear-gradient(90deg, #b4d8ff 0%, #9dc2ff 35%, #c5e4ff 100%);'>",
              "<span style='font-size: 28px; font-family: cursive; color: #000000;'>Categorical Variable V", var, "</span>",
              "</div>",
              "</div>"
            )),
            tableOutput(paste0("tablelam", var)),
            HTML("<br>"),
            HTML("<br>")
          ))
        }
        ui_elements
      })  
      
    }
    
    else if (choice == "Ising Model") {
      output$NuuTables <- renderUI({
        if (is.null(selected_binary_dependent_covariates)) return(NULL)
        
        ui_elements <- list()
        ui_elements <- append(ui_elements, list(
            HTML(paste0(
              "<div style='text-align: center; padding: 20px;'>",
              "<div style='border: 2px solid black; background-color: #FFAC1C; padding: 5px; display: inline-block; border-radius: 10px; box-shadow: 0px 0px #000000; background-image: linear-gradient(90deg, #b4d8ff 0%, #9dc2ff 35%, #c5e4ff 100%);'>",
              "<span style='font-size: 28px; font-family: cursive; color: #000000;'>Threshold Parameters", "</span>",
              "</div>",
              "</div>"
            )),
            tableOutput("tablenuu"),
            HTML("<br>"),
            HTML("<br>")
          ))
        
        ui_elements
      })
      
      output$IntPlot <- renderUI({
        if (is.null(selected_binary_dependent_covariates)) return(NULL)
        ui_elements <- list()
        ui_elements <- append(ui_elements, list(
          HTML(paste0(
            "<div style='text-align: center; padding: 20px;'>",
            "<div style='border: 2px solid black; background-color: #FFAC1C; padding: 5px; display: inline-block; border-radius: 10px; box-shadow: 0px 0px #000000; background-image: linear-gradient(90deg, #b4d8ff 0%, #9dc2ff 35%, #c5e4ff 100%);'>",
            "<span style='font-size: 28px; font-family: cursive; color: #000000;'>Interaction Parameters", "</span>",
            "</div>",
            "</div>"
          ))
        ))
        for(c in 1:C){
          ui_elements <- append(ui_elements, list(
            plotOutput(paste0("plot", c)),
            HTML("<br>"),
            HTML("<br>")
          ))
        }
        
        ui_elements
      })
      
    }
    
    
  })
  
  
  #--------------------------CONTINUOUS COVARIATES------------------------------
  output$tablemu <- renderText({
    if (input$choice == "Continuous Covariates Parameters") {
      mu <- data.frame(matrix(ncol = 0,nrow = dim(U)[2]))
      mu$Variable <- names(fitting$params$mu[,,1])
      for(c in 1:C){
        coln <- paste0("mu","_",c)
        mu[[coln]] = signif(fitting$params$mu[,,c], digits = 4)
      }
      colnames(mu) <- c("Variable", paste0("\u03BC",pedix[1:C]) ) 
      mu1 <- as.matrix(mu)
      colnames(mu1) <- NULL
      kable_mu <- mu1 %>% kbl(align = c("c","c","c","c")) %>%
        kable_styling(font_size = 20) %>% 
        column_spec(1, color = "black",
                    background = alpha("grey",0.2), bold = TRUE)  %>% 
        column_spec(2, color = "black",
                    background = my_colors[1]) %>%
        column_spec(3, color = "black",
                    background = my_colors[2]) 
      
      if(C == 3 | C == 4 | C == 5 | C == 6){
        kable_mu <- kable_mu %>% column_spec(4, color = "black",
                                             background = my_colors[3])
      }
      if(C == 4 | C == 5 | C == 6){
        kable_mu <- kable_mu %>% column_spec(5, color = "black",
                                             background = my_colors[4])
      }
      if(C == 5 | C == 6){
        kable_mu <- kable_mu %>% column_spec(6, color = "black",
                                             background = my_colors[5])
      }
      if(C == 6){
        kable_mu <- kable_mu %>% column_spec(7, color = "black",
                                             background = my_colors[6])
      }
      
      kable_mu <- kable_mu %>% column_spec(1:(C+1), border_right = "2px solid black", border_left = "2px solid black") %>%
        row_spec(seq(1,dim(U)[2]), extra_css = "border-bottom: 2px solid;") %>%
        column_spec(1, extra_css = "border-left: 2px solid;") %>%
        add_header_above(colnames(mu), extra_css = border_string) 
      return(kable_mu)
    }
  })
  
  output$tablesigma <- renderText({
    if (input$choice == "Continuous Covariates Parameters") {
      sigma <- data.frame(matrix(ncol = 0,nrow = dim(U)[2]))
      var_names <- names(fitting$params$mu[,,1])
      sigma$Variable <- var_names
      for(c in 1:C){
        sigma <- cbind(sigma,signif(fitting$params$sigma[,,c],digits = 4),
                       names(fitting$params$mu[,,c]))
      }
      sigma <- sigma[,-ncol(sigma)]
      sig_vec <- paste0("\u03A3",pedix[1:C])
      vec <- NULL
      for(i in 1:C){
        vec <- c(vec,"Variables", rep(sig_vec[i],dim(U)[2]) )
      }
      vec2 <- NULL
      for(i in 1:C){
        vec2 <- c(vec2," ", rep(sig_vec[i],dim(U)[2]) )
      }
      colnames(sigma) <- vec
      tt <- data.frame(matrix(  rep(c("Variables",
                                var_names),C) ,ncol = dim(U)[2]*C + C,
                                nrow = 1)  )
      colnames(tt) <- vec
      sigma2 <- rbind(tt,sigma)
      sigma1 <- as.matrix(sigma2)
      colnames(sigma2) <- NULL
      rownames(sigma2) <- NULL
      step_size <- dim(U)[2] - 1
      kable_sigma <- sigma2 %>% kbl(align = c("c","c","c","c")) %>%
        kable_styling(font_size = 20) %>% column_spec(seq(1,dim(U)[2]*C + C,dim(U)[2]+1), color = "black",
                                  background = alpha("grey",0.2), bold = TRUE)  %>%
        column_spec(2:(2+step_size), color = "black",
                    background = my_colors[1]) %>%
        column_spec((4+step_size):(4+2*step_size), color = "black",
                    background = my_colors[2])
      
      if(C == 3 | C == 4 | C == 5 | C == 6){
        kable_sigma <- kable_sigma %>% column_spec((6+2*step_size):(6+3*step_size),
                                      color = "black",
                                      background = my_colors[3])
      }
      if(C == 4 | C == 5 | C == 6){
        kable_sigma <- kable_sigma %>% column_spec((8+3*step_size):(8+4*step_size), color = "black",
                                                   background = my_colors[4])
      }
      if(C == 5 | C == 6){
        kable_sigma <- kable_sigma %>% column_spec((10+4*step_size):(10+5*step_size), color = "black",
                                                   background = my_colors[5])
      }
      if(C == 6){
        kable_sigma <- kable_sigma %>% column_spec((12+5*step_size):(12+6*step_size), color = "black",
                                                   background = my_colors[6])
      }
      
      kable_sigma <- kable_sigma %>% row_spec(1, color = "black",
                                              background = alpha("grey",0.2), bold = TRUE) %>%
        column_spec(seq(1,dim(U)[2]*C + C,dim(U)[2]+1), border_right = "2px solid black", border_left = "2px solid black") %>%
        column_spec(dim(U)[2]*C + C, border_right = "2px solid black") %>%
        row_spec(c(1,dim(U)[2]+1), extra_css = "border-bottom: 2px solid;") %>%
        add_header_above(vec2)
      return(kable_sigma)
    }
  })
  
  #--------------------------CATEGORICAL COVARIATES------------------------------
  output$tablelam1 <- renderText({
    if (input$choice == "Categorical Covariates Parameters") {
    V1 <- V[[1]]
    fir <- data.frame(matrix(ncol = 0,nrow = dim(V1)[2]))
    fir$Category <- sub(".*\\.", "", colnames(V1))
    fir1 <- data.frame(cbind(lambda_11 =  signif(fitting$params$lambda[[1]][,,1],digits = 4),
                            lambda_12 = signif(fitting$params$lambda[[1]][,,2],digits = 4)) )
    if(C == 3 | C == 4 | C == 5 | C == 6){
      fir1 <- cbind(fir1, lambda_13 = signif(fitting$params$lambda[[1]][,,3],digits = 4))
    }
    if(C == 4 | C == 5 | C == 6){
      fir1 <- cbind(fir1, lambda_14 = signif(fitting$params$lambda[[1]][,,4],digits = 4))
    }
    if(C == 5 | C == 6){
      fir1 <- cbind(fir1, lambda_15 = signif(fitting$params$lambda[[1]][,,5],digits = 4))
    }
    if(C == 6){
      fir1 <- cbind(fir1, lambda_16 = signif(fitting$params$lambda[[1]][,,6],digits = 4))
    }
    
    fir <- cbind(fir,fir1)
    colnames(fir) <- c("Category", paste0("\u03BB",pedix[1], pedix[1:C]) )
    fir1 <- as.matrix(fir)
    colnames(fir1) <- NULL
    rownames(fir1) <- NULL
    kable_fir1 <- fir1 %>% kbl(align = c("c","c","c","c","c")) %>%
      kable_styling(font_size = 20) %>% 
     column_spec(1, color = "black",
                      background = alpha("grey",0.2)) %>%
      column_spec(2, color = "black",
         background = my_colors[1]) %>%
      column_spec(3, color = "black",
                  background = my_colors[2])
      
      if(C == 3 | C == 4  | C == 5 | C == 6){
        kable_fir1 <- kable_fir1  %>% column_spec(4, color = "black",
                                             background = my_colors[3])
      }
       if(C == 4| C == 5 | C == 6){
         kable_fir1 <- kable_fir1  %>% column_spec(5, color = "black",
                                           background = my_colors[4])
       }
      if(C == 5 | C == 6){
      kable_fir1 <- kable_fir1  %>% column_spec(6, color = "black",
                                                background = my_colors[5])
      }
    if(C == 6){
      kable_fir1 <- kable_fir1  %>% column_spec(7, color = "black",
                                                background = my_colors[6])
    }
    
    kable_fir1 <- kable_fir1 %>%
      column_spec(seq(1,C+1), border_right = "2px solid black", border_left = "2px solid black") %>%
      row_spec(seq(1,dim(V1)[2]), extra_css = "border-bottom: 2px solid;") %>%
      add_header_above(colnames(fir), extra_css = border_string) 
    
    return(kable_fir1)
    }
  })
  
  
  output$tablelam2 <- renderText({
    if (input$choice == "Categorical Covariates Parameters") {
      V2 <- V[[2]]
      sec <- data.frame(matrix(ncol = 0,nrow = dim(V2)[2]))
      sec$Category <- sub(".*\\.", "", colnames(V2))
      sec1 <- data.frame(cbind(lambda_21 = 
                                signif(fitting$params$lambda[[2]][,,1], digits = 4),
                              lambda_22 = signif(fitting$params$lambda[[2]][,,2], digits = 4) ) )
      if(C == 3 | C == 4 | C == 5 | C == 6){
        sec1 <- cbind(sec1, lambda_13 = signif(fitting$params$lambda[[2]][,,3],digits = 4) )
      }
      if(C == 4 | C == 5 | C == 6){
        sec1 <- cbind(sec1, lambda_14 = signif(fitting$params$lambda[[2]][,,4],digits = 4) )
      }
      if(C == 5 | C == 6){
        sec1 <- cbind(sec1, lambda_15 = signif(fitting$params$lambda[[2]][,,5],digits = 4) )
      }
      if(C == 6){
        sec1 <- cbind(sec1, lambda_16 = signif(fitting$params$lambda[[2]][,,6],digits = 4) )
      }
      sec <- cbind(sec,sec1)
      colnames(sec) <- c("Category", paste0("\u03BB",pedix[2], pedix[1:C]) )
      sec1 <- as.matrix(sec)
      colnames(sec1) <- NULL
      rownames(sec1) <- NULL
      kable_sec1 <- sec1 %>% kbl(align = c("c","c","c","c")) %>%
        kable_styling(font_size = 20) %>%
        column_spec(1, color = "black",
                    background = alpha("grey",0.2)) %>%
        column_spec(2, color = "black",
                background = my_colors[1]) %>%
        column_spec(3, color = "black",
                    background = my_colors[2]) 
      
        if(C == 3 | C == 4 | C == 5 | C == 6){
          kable_sec1 <- kable_sec1 %>% column_spec(4, color = "black",
                                       background = my_colors[3])
        }
      if(C == 4 | C == 5 | C == 6){
        kable_sec1 <- kable_sec1 %>% column_spec(5, color = "black",
                                     background = my_colors[4])
      }
      if(C == 5 | C == 6){
        kable_sec1 <- kable_sec1 %>% column_spec(6, color = "black",
                                                 background = my_colors[5])
      }
      if(C == 6){
        kable_sec1 <- kable_sec1 %>% column_spec(7, color = "black",
                                                 background = my_colors[6])
      }
      kable_sec1 <- kable_sec1 %>%
        column_spec(seq(1,C+1), border_right = "2px solid black", border_left = "2px solid black") %>%
        row_spec(seq(1,dim(V2)[2]), extra_css = "border-bottom: 2px solid;") %>%
        add_header_above(colnames(sec), extra_css = border_string) 
      return(kable_sec1)
    }
  })
  
  
  output$tablelam3 <- renderText({
    if (input$choice == "Categorical Covariates Parameters") {
      V3 <- V[[3]]
      thi <- data.frame(matrix(ncol = 0,nrow = dim(V3)[2]))
      thi$Category <- 1:dim(V3)[2]
      thi1 <- data.frame(cbind(lambda_31 = 
                                signif(fitting$params$lambda[[3]][,,1], digits = 4),
                              lambda_32 = signif(fitting$params$lambda[[3]][,,2] , digits = 4) ))
      if(C == 3 | C == 4 | C == 5 | C == 6){
        thi1 <- cbind(thi1, lambda_33 = signif(fitting$params$lambda[[3]][,,3],digits = 4))
      }
      if(C == 4 | C == 5 | C == 6){
        thi1 <- cbind(thi1, lambda_34 = signif(fitting$params$lambda[[3]][,,4],digits = 4))
      }
      if(C == 5 | C == 6){
        thi1 <- cbind(thi1, lambda_35 = signif(fitting$params$lambda[[3]][,,5],digits = 4))
      }
      if(C == 6){
        thi1 <- cbind(thi1, lambda_36 = signif(fitting$params$lambda[[3]][,,6],digits = 4))
      }
      thi <- cbind(thi,thi1)
      colnames(thi) <- c("Category", paste0("\u03BB",pedix[3], pedix[1:C]) )
      thi1 <- as.matrix(thi)
      colnames(thi1) <- NULL
      kable_thi1 <- thi1 %>% kbl(align = c("c","c","c","c")) %>%
        kable_styling(font_size = 20) %>% 
        column_spec(1, color = "black",
                    background = alpha("grey",0.2)) %>%
        column_spec(2, color = "black",
                  background = my_colors[1]) %>%
        column_spec(3, color = "black",
                    background = my_colors[2])
      
        if(C == 3 | C == 4 | C == 5 | C == 6){
          kable_thi1 <- kable_thi1 %>% column_spec(4, color = "black",
                                       background = my_colors[3])
        }
      if(C == 4 | C == 5 | C == 6){
        kable_thi1 <- kable_thi1 %>% column_spec(5, color = "black",
                                     background = my_colors[4])
      }
      if(C == 5 | C == 6){
        kable_thi1 <- kable_thi1 %>% column_spec(6, color = "black",
                                                 background = my_colors[5])
      }
      if(C == 6){
        kable_thi1 <- kable_thi1 %>% column_spec(7, color = "black",
                                                 background = my_colors[6])
      }
      
      kable_thi1 <- kable_thi1 %>%
        column_spec(seq(1,C+1), border_right = "2px solid black", border_left = "2px solid black") %>%
        row_spec(seq(1,dim(V3)[2]), extra_css = "border-bottom: 2px solid;") %>%
        add_header_above(colnames(thi), extra_css = border_string) 
      return(kable_thi1)
    }
  })
  
  
 
  
  #--------------------------FIXED EFFECTS--------------------------------------
  output$tablefix1 <- renderText({
    if (input$choice == "Fixed Effects") {
      fix1 <- as.data.frame( signif(summary(fitting$params$models[[1]])$coefficients,digits = 3) )
      fixa <- as.matrix(fix1)
      colnames(fixa) <- NULL
      fixa %>% kbl(align = c("c","c","c","c","c","c")) %>%
        kable_styling(font_size = 15) %>% 
        row_spec(which(fix1[,4] < 0.05),
                 color = "black",
                 background = "rgba(0, 0, 255, 0.5)") %>%
        row_spec(which(fix1[,4] >= 0.05),
                 color = "black",
                 background = alpha("grey", 0.3)) %>%
        column_spec(1, color = "black", bold = TRUE) %>%
        column_spec(seq(1,5), border_right = "2px solid black", border_left = "2px solid black") %>%
        row_spec(seq(1, nrow(summary(fitting$params$models[[1]])$coefficients)), 
                 extra_css = "border-bottom: 2px solid;") %>%
        add_header_above(c("Variables",colnames(fix1)), extra_css = border_string) 
    }
  })
  
  
  output$tablefix2 <- renderText({
    if (input$choice == "Fixed Effects") {
      fix1 <- as.data.frame( signif(summary(fitting$params$models[[2]])$coefficients,digits = 3) )
      fixa <- as.matrix(fix1)
      colnames(fixa) <- NULL
      fixa %>% kbl(align = c("c","c","c","c","c","c")) %>%
        kable_styling(font_size = 15) %>% 
        row_spec(which(fix1[,4] < 0.05),
                 color = "black",
                 background = "rgba(0, 255, 0, 0.4)") %>%
        row_spec(which(fix1[,4] >= 0.05),
                 color = "black",
                 background = alpha("grey", 0.3)) %>%
        column_spec(1, color = "black", bold = TRUE) %>%
        column_spec(seq(1,5), border_right = "2px solid black", border_left = "2px solid black") %>%
        row_spec(seq(1, nrow(summary(fitting$params$models[[1]])$coefficients)), extra_css = "border-bottom: 2px solid;")  %>%
        add_header_above(c("Variables",colnames(fix1)), extra_css = border_string) 
    }
  })
  
  
  output$tablefix3 <- renderText({
    if (input$choice == "Fixed Effects") {
      fix1 <- as.data.frame( signif(summary(fitting$params$models[[3]])$coefficients,digits = 3) )
      fixa <- as.matrix(fix1)
      colnames(fixa) <- NULL
      fixa %>% kbl(align = c("c","c","c","c","c","c")) %>%
        kable_styling(font_size = 15) %>% 
        row_spec(which(fix1[,4] < 0.05),
                 color = "black",
                 background = my_colors[3]) %>%
        row_spec(which(fix1[,4] >= 0.05),
                 color = "black",
                 background = alpha("grey", 0.3)) %>%
        column_spec(1, color = "black", bold = TRUE) %>%
        column_spec(seq(1,5), border_right = "2px solid black", border_left = "2px solid black") %>%
        row_spec(seq(1, nrow(summary(fitting$params$models[[1]])$coefficients)), extra_css = "border-bottom: 2px solid;") %>%
        add_header_above(c("Variables",colnames(fix1)), extra_css = border_string) 
    }
  })
  
  
  output$tablefix4 <- renderText({
    if (input$choice == "Fixed Effects") {
      fix1 <- as.data.frame( signif(summary(fitting$params$models[[4]])$coefficients,digits = 3) )
      fixa <- as.matrix(fix1)
      colnames(fixa) <- NULL
      fixa %>% kbl(align = c("c","c","c","c","c","c")) %>%
        kable_styling(font_size = 15) %>% 
        row_spec(which(fix1[,4] < 0.05),
                 color = "black",
                 background = my_colors[4]) %>%
        row_spec(which(fix1[,4] >= 0.05),
                 color = "black",
                 background = alpha("grey", 0.3)) %>%
        column_spec(1, color = "black", bold = TRUE) %>%
        column_spec(seq(1,5), border_right = "2px solid black", border_left = "2px solid black") %>%
        row_spec(seq(1,dim(U)[2] + dim(D)[2] +
                       ifelse(length(V) >= 1, dim(V[[1]])[2] - 1, 0) +
                       ifelse(length(V) >= 2, dim(V[[2]])[2] - 1, 0) +
                       ifelse(length(V) >= 3, dim(V[[3]])[2] - 1, 0)), extra_css = "border-bottom: 2px solid;") %>%
        add_header_above(c("Variables",colnames(fix1)), extra_css = border_string) 
    }
  })
  
  output$tablefix5 <- renderText({
    if (input$choice == "Fixed Effects") {
      fix1 <- as.data.frame( signif(summary(fitting$params$models[[5]])$coefficients,digits = 3) )
      fixa <- as.matrix(fix1)
      colnames(fixa) <- NULL
      fixa %>% kbl(align = c("c","c","c","c","c","c")) %>%
        kable_styling(font_size = 15) %>% 
        row_spec(which(fix1[,4] < 0.05),
                 color = "black",
                 background = my_colors[5]) %>%
        row_spec(which(fix1[,4] >= 0.05),
                 color = "black",
                 background = alpha("grey", 0.3)) %>%
        column_spec(1, color = "black", bold = TRUE) %>%
        column_spec(seq(1,5), border_right = "2px solid black", border_left = "2px solid black") %>%
        row_spec(seq(1,dim(U)[2] + dim(D)[2] +
                       ifelse(length(V) >= 1, dim(V[[1]])[2] - 1, 0) +
                       ifelse(length(V) >= 2, dim(V[[2]])[2] - 1, 0) +
                       ifelse(length(V) >= 3, dim(V[[3]])[2] - 1, 0)), extra_css = "border-bottom: 2px solid;") %>%
        add_header_above(c("Variables",colnames(fix1)), extra_css = border_string) 
    }
  })
  
  output$tablefix6 <- renderText({
    if (input$choice == "Fixed Effects") {
      fix1 <- as.data.frame( signif(summary(fitting$params$models[[6]])$coefficients,digits = 3) )
      fixa <- as.matrix(fix1)
      colnames(fixa) <- NULL
      fixa %>% kbl(align = c("c","c","c","c","c","c")) %>%
        kable_styling(font_size = 15) %>% 
        row_spec(which(fix1[,4] < 0.05),
                 color = "black",
                 background = my_colors[6]) %>%
        row_spec(which(fix1[,4] >= 0.05),
                 color = "black",
                 background = alpha("grey", 0.3)) %>%
        column_spec(1, color = "black", bold = TRUE) %>%
        column_spec(seq(1,5), border_right = "2px solid black", border_left = "2px solid black") %>%
        row_spec(seq(1,dim(U)[2] + dim(D)[2] +
                       ifelse(length(V) >= 1, dim(V[[1]])[2] - 1, 0) +
                       ifelse(length(V) >= 2, dim(V[[2]])[2] - 1, 0) +
                       ifelse(length(V) >= 3, dim(V[[3]])[2] - 1, 0)), extra_css = "border-bottom: 2px solid;") %>%
        add_header_above(c("Variables",colnames(fix1)), extra_css = border_string) 
    }
  })
  
  
  #----------------------------BIC COMPARISON-----------------------------------
  output$biccomp <- renderText({
      border_string <- "font-family: cursive; border-bottom: 2px solid;border-top: 2px solid;border-left: 2px solid;border-right: 2px solid;"
      ff = as.matrix(cbind(Clusters = paste("C =", C_vec), BIC = round(BIC_vec,2)))
      cn <- colnames(ff)
      idx_best <- which.min(BIC_vec)
      colnames(ff) <- NULL
      ff %>% kbl(align = c("c", "c")) %>%
        kable_styling(font_size = 20)  %>%
        row_spec(idx_best, color = "black",
                 background = "rgba(0, 255, 0, 0.6)", bold = F) %>%
        row_spec(seq(1,length(BIC_vec))) %>%
        column_spec(1:2) %>%
        add_header_above(cn, extra_css = border_string, color = "white", 
                         background = alpha("#42B3FF",0.2)) 
  })
  

  
  
  
  #--------------------------ISING MODEL----------------------------------------
  output$tablenuu <- renderText({
    if (input$choice == "Ising Model") {
      
    nuu <- data.frame(matrix(ncol = 0,nrow = dim(D)[2]  ))
    nuu$Variable <- paste0("D",1:dim(D)[2])
    for(c in 1:C){
      coln <- paste0("nu","\u208",c)
      nuu[[coln]] = signif(fitting$params$thres[,,c], digits = 3)
    }
    colnames(nuu) <- c("Variable", paste0("\u03BD",pedix[1:C]) )
    nuu1 <- as.matrix(nuu)
    colnames(nuu1) <- NULL
    kable_nuu <- nuu1 %>% kbl(align = c("c","c","c")) %>%
      kable_styling(font_size = 20) %>% column_spec(1, color = "black",
                                                  background = alpha("grey",0.2)) %>%
      column_spec(2, color = "black",
                background = my_colors[1]) %>%
      column_spec(3, color = "black",
                background = my_colors[2]) 
      
      if(C == 3 | C == 4 | C == 5 | C == 6){
        kable_nuu <- kable_nuu %>% column_spec(4, color = "black",
                                                 background = my_colors[3])
      }
    if(C == 4 | C == 5 | C == 6){
      kable_nuu <- kable_nuu %>% column_spec(5, color = "black",
                                               background = my_colors[4])
    }
    if(C == 5 | C == 6){
      kable_nuu <- kable_nuu %>% column_spec(6, color = "black",
                                             background = my_colors[5])
    }
    if(C == 6){
      kable_nuu <- kable_nuu %>% column_spec(7, color = "black",
                                             background = my_colors[6])
    }
    
    kable_nuu <- kable_nuu %>%
      column_spec(1:(C+1), border_right = "2px solid black", border_left = "2px solid black") %>%
      row_spec(seq(1,dim(D)[2]), extra_css = "border-bottom: 2px solid;") %>%
      add_header_above(colnames(nuu), extra_css = border_string) 
    }
  })
  
  output$plot1 <- renderPlot({
    if (input$choice == "Ising Model") {
      lab <- paste0("D",1:dim(D)[2])
      qgraph(fitting$params$int[,,1], fade = FALSE, 
             edge.labels = round(fitting$params$int[,,1], 3), 
             edge.label.cex = 1.5,labels = lab, edge.label.color = "black",
             color = alpha("blue", 0.5),edge.width = 0.7, 
             negCol = "red",posCol = "forestgreen",title = "First cluster Fitted network",
             title.cex = 2.5,edge.label.margin = 0.001,label.font = 2,edge.label.font = 2) 
    }
  })
  output$plot2 <- renderPlot({
    if (input$choice == "Ising Model") {
      lab <- paste0("D",1:dim(D)[2])
      qgraph(fitting$params$int[,,2], fade = FALSE, edge.labels = round(fitting$params$int[,,2], 3), 
             edge.label.cex = 1.5,labels = lab, edge.label.color = "black",color = alpha("green", 0.5),
             edge.width = 0.7, negCol = "red",posCol = "forestgreen",weighted = T,
             title = "Second cluster Fitted network",title.cex =2.5,edge.label.margin = 0.001,
             label.font = 2,edge.label.font = 2)
    }
  })
  output$plot3 <- renderPlot({
    if (input$choice == "Ising Model") {
      lab <- paste0("D",1:dim(D)[2])
      qgraph(fitting$params$int[,,3], fade = FALSE, 
             edge.labels = round(fitting$params$int[,,3], 3), 
             edge.label.cex = 1.5,labels = lab, edge.label.color = "black",
             color = alpha("deeppink", 0.5),edge.width = 0.7
             , negCol = "red",posCol = "forestgreen",weighted = T,title = "Third cluster Fitted network",
             title.cex = 2.5,edge.label.margin = 0.001,label.font = 2,edge.label.font = 2)
    }
  })
  
  output$plot4 <- renderPlot({
    if (input$choice == "Ising Model") {
      lab <- paste0("D",1:dim(D)[2])
      qgraph(fitting$params$int[,,4], fade = FALSE,
             edge.labels = round(fitting$params$int[,,4], 3),
             edge.label.cex = 1.5,labels = lab, edge.label.color = "black",
             color = alpha("orange", 0.5),edge.width = 0.7
             , negCol = "red",posCol = "forestgreen",weighted = T,title = "Fourth cluster Fitted network",
             title.cex = 2.5,edge.label.margin = 0.001,label.font = 2,edge.label.font = 2)
    }
  })
  
  output$plot5 <- renderPlot({
    if (input$choice == "Ising Model") {
      lab <- paste0("D",1:dim(D)[2])
      qgraph(fitting$params$int[,,5], fade = FALSE,
             edge.labels = round(fitting$params$int[,,5], 3),
             edge.label.cex = 1.5,labels = lab, edge.label.color = "black",
             color = alpha("red", 0.5),edge.width = 0.7
             , negCol = "red",posCol = "forestgreen",weighted = T,title = "Fifth cluster Fitted network",
             title.cex = 2.5,edge.label.margin = 0.001,label.font = 2,edge.label.font = 2)
    }
  })
  
  output$plot6 <- renderPlot({
    if (input$choice == "Ising Model") {
      lab <- paste0("D",1:dim(D)[2])
      qgraph(fitting$params$int[,,6], fade = FALSE,
             edge.labels = round(fitting$params$int[,,5], 3),
             edge.label.cex = 1.5,labels = lab, edge.label.color = "black",
             color = alpha("purple", 0.5),edge.width = 0.7
             , negCol = "red",posCol = "forestgreen",weighted = T,title = "Sixth cluster Fitted network",
             title.cex = 2.5,edge.label.margin = 0.001,label.font = 2,edge.label.font = 2)
    }
  })
  
  #--------------------------RANDOM EFFECTS-------------------------------------
  output$tableranef <- renderText({
    if (input$choice == "Random Effects") {
      
      ran_eff <- data.frame(matrix(ncol = 0,nrow = length( levels(as.factor(mydata$level)) )  ))
      ran_eff$Levels <- levels(as.factor(mydata$level))
      for(c in 1:C){
        coln <- paste0("b","\u208",c)
        ran_eff[[coln]] = signif(ranef(fitting$params$models[[c]])$level[,1], digits = 3)
      }
      colnames(ran_eff) <- c("Level",paste0("b",pedix[1:C]))
      ran_eff1 <- as.matrix(ran_eff)
      colnames(ran_eff1) <- NULL
      kable_ranef <- ran_eff1 %>% kbl(align = c("c","c","c")) %>%
        kable_styling(font_size = 20) %>% column_spec(1, color = "black",
                                                      background = alpha("grey",0.2)) %>%
        column_spec(2, color = "black",
                    background = my_colors[1]) %>%
        column_spec(3, color = "black",
                    background = my_colors[2]) 
        
        if(C == 3 | C == 4 | C == 5 | C == 6){
          kable_ranef <- kable_ranef %>% column_spec(4, color = "black",
                                       background = my_colors[3])
        }
      if(C == 4 | C == 5 | C == 6){
        kable_ranef <- kable_ranef %>% column_spec(5, color = "black",
                                     background = my_colors[4])
      }
      if(C == 5 | C == 6){
        kable_ranef <- kable_ranef %>% column_spec(6, color = "black",
                                                   background = my_colors[5])
      }
      if(C == 6){
        kable_ranef <- kable_ranef %>% column_spec(7, color = "black",
                                                   background = my_colors[6])
      }
      kable_ranef <- kable_ranef %>%
        column_spec(1:(C+1), border_right = "2px solid black", border_left = "2px solid black") %>%
        row_spec(seq(1,length( levels(as.factor(mydata$level)) )), 
                 extra_css = "border-bottom: 2px solid;") %>%
        add_header_above(colnames(ran_eff), extra_css = border_string) 
      
      return(kable_ranef)
    }
  })
  
  output$plotre1 <- renderPlot({
    if (input$choice == "Random Effects") {
      re_sim <- REsim(fitting$params$models[[1]])
      re_sim$CI_lower <- re_sim$mean - (1.96 * re_sim$sd)
      re_sim$CI_upper <- re_sim$mean + (1.96 * re_sim$sd)
      re_sim$contains_zero <- with(re_sim, CI_lower <= 0 & CI_upper >= 0)
      re_sim$groupFctr <- " "
      re_sim <- re_sim %>%
        arrange(median)
      x_colors <- ifelse(re_sim$contains_zero == TRUE, "gray30", "blue")
      plotREsim(re_sim,labs = TRUE) +
        geom_hline(yintercept = 0, color = I("black"), size = 1.2) +
        geom_errorbar(aes(color = contains_zero), width = 0, size = 0.5) +
        geom_point(aes(color = contains_zero), size = 3) +
        scale_color_manual(values = c("blue", "gray50"), guide = "none") +
        theme_minimal() +
        theme(legend.position = "none",axis.text.x = element_text(size = 16,face = "bold",
                                                                  vjust = 3,angle = 0, color = x_colors),
              plot.title = element_text(color = "blue",size = 20, vjust = -2),
              axis.title = element_text(size = 18),
              axis.text.y = element_text(size = 16,face = "bold", hjust = 1.5))
      
    }
  })
  output$plotre2 <- renderPlot({
    if (input$choice == "Random Effects") {
      re_sim <- REsim(fitting$params$models[[2]])
      re_sim$CI_lower <- re_sim$mean - (1.96 * re_sim$sd)
      re_sim$CI_upper <- re_sim$mean + (1.96 * re_sim$sd)
      re_sim$contains_zero <- with(re_sim, CI_lower <= 0 & CI_upper >= 0)
      re_sim$groupFctr <- " "
      re_sim <- re_sim %>%
        arrange(median)
      x_colors <- ifelse(re_sim$contains_zero == TRUE, "grey30", "green")
      plotREsim(re_sim,labs = TRUE) +
        geom_hline(yintercept = 0, color = I("black"), size = 1.2) +
        geom_errorbar(aes(color = contains_zero), width = 0, size = 0.8) +
        geom_point(aes(color = contains_zero), size = 3) +
        scale_color_manual(values = c("green", "gray50"), guide = "none") +
        theme_minimal() +
        theme(legend.position = "none",axis.text.x = element_text(size = 16,face = "bold",
                                                                  vjust = 3,angle = 0, color = x_colors),
              plot.title = element_text(color = "green",size = 20, vjust = -2),
              axis.title = element_text(size = 18),
              axis.text.y = element_text(size = 16,face = "bold", hjust = 1.5))
      
    }
  })
  output$plotre3 <- renderPlot({
    if (input$choice == "Random Effects") {
      re_sim <- REsim(fitting$params$models[[3]])
      re_sim$CI_lower <- re_sim$mean - (1.96 * re_sim$sd)
      re_sim$CI_upper <- re_sim$mean + (1.96 * re_sim$sd)
      re_sim$contains_zero <- with(re_sim, CI_lower <= 0 & CI_upper >= 0)
      re_sim$groupFctr <- " "
      re_sim <- re_sim %>%
        arrange(median)
      x_colors <- ifelse(re_sim$contains_zero == TRUE, "grey30", "deeppink")
      plotREsim(re_sim,labs = TRUE) +
        geom_hline(yintercept = 0, color = I("black"), size = 1.2) +
        geom_errorbar(aes(color = contains_zero), width = 0, size = 0.8) +
        geom_point(aes(color = contains_zero), size = 3) +
        scale_color_manual(values = c("deeppink", "gray50"), guide = "none") +
        theme_minimal() +
        theme(legend.position = "none",axis.text.x = element_text(size = 16,face = "bold",
                                                                  vjust = 3,angle = 0, color = x_colors),
              plot.title = element_text(color = "deeppink",size = 20, vjust = -2),
              axis.title = element_text(size = 18),
              axis.text.y = element_text(size = 16,face = "bold", hjust = 1.5))
    }
  })
  
  output$plotre4 <- renderPlot({
    if (input$choice == "Random Effects") {
      re_sim <- REsim(fitting$params$models[[4]])
      re_sim$CI_lower <- re_sim$mean - (1.96 * re_sim$sd)
      re_sim$CI_upper <- re_sim$mean + (1.96 * re_sim$sd)
      re_sim$contains_zero <- with(re_sim, CI_lower <= 0 & CI_upper >= 0)
      re_sim$groupFctr <- " "
      re_sim <- re_sim %>%
        arrange(median)
      x_colors <- ifelse(re_sim$contains_zero == TRUE, "grey30", "orange")
      plotREsim(re_sim,labs = TRUE) +
        geom_hline(yintercept = 0, color = I("black"), size = 1.2) +
        geom_errorbar(aes(color = contains_zero), width = 0, size = 0.8) +
        geom_point(aes(color = contains_zero), size = 3) +
        scale_color_manual(values = c("orange", "gray50"), guide = "none") +
        theme_minimal() +
        theme(legend.position = "none",axis.text.x = element_text(size = 16,face = "bold",
                                                                  vjust = 3,angle = 0, color = x_colors),
              plot.title = element_text(color = "orange",size = 20, vjust = -2),
              axis.title = element_text(size = 18),
              axis.text.y = element_text(size = 16,face = "bold", hjust = 1.5))
      
    }
  })
  
  output$plotre5 <- renderPlot({
    if (input$choice == "Random Effects") {
      re_sim <- REsim(fitting$params$models[[5]])
      re_sim$CI_lower <- re_sim$mean - (1.96 * re_sim$sd)
      re_sim$CI_upper <- re_sim$mean + (1.96 * re_sim$sd)
      re_sim$contains_zero <- with(re_sim, CI_lower <= 0 & CI_upper >= 0)
      re_sim$groupFctr <- " "
      re_sim <- re_sim %>%
        arrange(median)
      x_colors <- ifelse(re_sim$contains_zero == TRUE, "grey30", "red")
      plotREsim(re_sim,labs = TRUE) +
        geom_hline(yintercept = 0, color = I("black"), size = 1.2) +
        geom_errorbar(aes(color = contains_zero), width = 0, size = 0.8) +
        geom_point(aes(color = contains_zero), size = 3) +
        scale_color_manual(values = c("red", "gray50"), guide = "none") +
        theme_minimal() +
        theme(legend.position = "none",axis.text.x = element_text(size = 16,face = "bold",
                                                                  vjust = 3,angle = 0, color = x_colors),
              plot.title = element_text(color = "red",size = 20, vjust = -2),
              axis.title = element_text(size = 18),
              axis.text.y = element_text(size = 16,face = "bold", hjust = 1.5))
      
    }
  })
  
  output$plotre6 <- renderPlot({
    if (input$choice == "Random Effects") {
      re_sim <- REsim(fitting$params$models[[6]])
      re_sim$CI_lower <- re_sim$mean - (1.96 * re_sim$sd)
      re_sim$CI_upper <- re_sim$mean + (1.96 * re_sim$sd)
      re_sim$contains_zero <- with(re_sim, CI_lower <= 0 & CI_upper >= 0)
      re_sim$groupFctr <- " "
      re_sim <- re_sim %>%
        arrange(median)
      x_colors <- ifelse(re_sim$contains_zero == TRUE, "grey30", "purple")
      plotREsim(re_sim,labs = TRUE) +
        geom_hline(yintercept = 0, color = I("black"), size = 1.2) +
        geom_errorbar(aes(color = contains_zero), width = 0, size = 0.8) +
        geom_point(aes(color = contains_zero), size = 3) +
        scale_color_manual(values = c("purple", "gray50"), guide = "none") +
        theme_minimal() +
        theme(legend.position = "none",axis.text.x = element_text(size = 16,face = "bold",
                                                                  vjust = 3,angle = 0, color = x_colors),
              plot.title = element_text(color = "purple",size = 20, vjust = -2),
              axis.title = element_text(size = 18),
              axis.text.y = element_text(size = 16,face = "bold", hjust = 1.5))
      
    }
  })

}
  

# Run the application
shinyApp(ui = ui, server = server)
