library(shiny)
library(dplyr)
library(kableExtra)
library(gridExtra)
library(ggplot2)


pedix <- c("\u2081","\u2082","\u2083","\u2084","\u2085","\u2086","\u2087","\u2088","\u2089")

b_rgb <- col2rgb("blue")
dp_rgb <- col2rgb("deeppink")
or_rgb <- col2rgb("orange")
g_rgb <- col2rgb("green")

my_colors <- c("rgba(0, 0, 255, 0.6)", "rgba(0, 255, 0, 0.6)",
               "rgba(255, 20, 147, 0.6)", "rgba(255, 165, 0, 0.6)")

#my_colors <- c(alpha("blue",0.4),my_colors[2],my_colors[3],
               #alpha("orange",0.4))
border_string <- "border-bottom: 2px solid;border-top: 2px solid;border-left: 2px solid;border-right: 2px solid;"

ui <- fluidPage(
  selectInput("choice", "Select a model parameter:", choices = c("Continuous Covariates Parameters",
                                                                 "Categorical Covariates Parameters",
                                                                 "Random Effects","Fixed Effects",
                                                                 "Ising Model")),
    tableOutput("outputTable"),
  conditionalPanel(
    condition = "input.choice == 'Ising Model'",
    fluidRow(
      column(4, plotOutput("plot1")),  # Adjust the width of the first plot
      column(4, plotOutput("plot2")),  # Adjust the width of the second plot
      column(4, plotOutput("plot3"))   # Adjust the width of the third plot
    )
    
  ),
  conditionalPanel(
    condition = "input.choice == 'Random Effects'",
    fluidRow(
      column(4, plotOutput("plotre1")),  # Adjust the width of the first plot
      column(4, plotOutput("plotre2")),  # Adjust the width of the second plot
      column(4, plotOutput("plotre3"))   # Adjust the width of the third plot
    )
    
  ),
  conditionalPanel(
    condition = "input.choice == 'Categorical Covariates Parameters'",
    uiOutput("conditionalUI")
  ),
  
  conditionalPanel(
    condition = "input.choice == 'Continuous Covariates Parameters'",
    HTML(paste0("<span style='font-size: 28px; font-weight: bold; font-family: Arial;'>Note:</span>","<p><span style='font-size: 20px; font-weight: normal; font-family: Arial;'>The continuous variables are modeled, in each cluster, with a multivariate normal distribution with mean <span style='vertical-align:sub;'>μ<sub>c</sub></span> and covariance matrix <span style='vertical-align:sub;'>Σ<sub>c</sub></span></span></p>")),
    HTML(paste0(
      "<div style='text-align: center; padding: 20px;'>",
      "<div style='border: 2px solid black; background-color: whitesmoke; padding: 10px; display: inline-block; border-radius: 10px;'>",
      "<span style='font-size: 28px; font-weight: bold; font-family: 'Arial', sans-serif; color: #0077B6;'>Continuous covariates mean", "</span>",
      "</div>",
      "</div>"
    )),
    tableOutput("tablemu"),
    HTML("<br>"),
    HTML(paste0(
      "<div style='text-align: center; padding: 20px;'>",
      "<div style='border: 2px solid black; background-color: whitesmoke; padding: 10px; display: inline-block; border-radius: 10px;'>",
      "<span style='font-size: 28px; font-weight: bold; font-family: 'Arial', sans-serif; color: #0077B6;'>Continuous covariates covariance matrix", "</span>",
      "</div>",
      "</div>"
    )),
    tableOutput("tablesigma"),
    HTML("<br>"),
    HTML("<br>"),
    plotOutput("plotclust")
  ),
  conditionalPanel(
    condition = "input.choice == 'Fixed Effects'",
    fluidRow(
      column(4, style = "padding-right: 10px;",tableOutput("tablefix1")),
    column(4, style = "padding-right: 10px;",tableOutput("tablefix2")),
    column(4, style = "padding-right: 10px;",tableOutput("tablefix3"))
    )
    )

  
)

server <- function(input, output) {
  output$outputTable <- renderText({
    choice <- input$choice
    if (choice == "Random Effects") {
      
      ran_eff <- data.frame(matrix(ncol = 0,nrow = length( levels(as.factor(data$level)) )  ))
      ran_eff$Levels <- levels(as.factor(data$level))
      for(c in 1:C){
        coln <- paste0("b","\u208",c)
        ran_eff[[coln]] = signif(ranef(fitting$params$models[[c]])$level[,1], digits = 3)
      }
      colnames(ran_eff) <- c("Level",paste0("b",pedix[1:C]))
      ran_eff1 <- as.matrix(ran_eff)
      colnames(ran_eff1) <- NULL
      ran_eff1 %>% kbl(align = c("c","c","c")) %>%
        kable_styling(font_size = 20) %>% column_spec(1, color = "black",
                                        background = alpha("grey",0.2)) %>%
        column_spec(2, color = "black",
                                        background = my_colors[1]) %>%
        column_spec(3, color = "black",
                    background = my_colors[2]) %>%
        column_spec(4, color = "black",
                    background = my_colors[3]) %>%
        column_spec(1:(C+1), border_right = "2px solid black", border_left = "2px solid black") %>%
        row_spec(seq(1,length( levels(as.factor(data$level)) )), 
                 extra_css = "border-bottom: 2px solid;") %>%
        add_header_above(colnames(ran_eff), extra_css = border_string) 
      
    }
    
    
    else if (choice == "Ising Model") {
      
      nuu <- data.frame(matrix(ncol = 0,nrow = dim(D)[2]  ))
      nuu$Variable <- paste0("D",1:dim(D)[2])
      for(c in 1:C){
        coln <- paste0("nu","\u208",c)
        nuu[[coln]] = signif(fitting$params$thres[,,c], digits = 3)
      }
      colnames(nuu) <- c("Variable", paste0("\u03BD",pedix[1:C]) )
      nuu1 <- as.matrix(nuu)
      colnames(nuu1) <- NULL
      nuu1 %>% kbl(align = c("c","c","c")) %>%
        kable_styling(font_size = 20) %>% column_spec(1, color = "black",
                                                      background = alpha("grey",0.2)) %>%
        column_spec(2, color = "black",
                    background = my_colors[1]) %>%
        column_spec(3, color = "black",
                    background = my_colors[2]) %>%
        column_spec(4, color = "black",
                    background = my_colors[3]) %>%
        column_spec(1:(C+1), border_right = "2px solid black", border_left = "2px solid black") %>%
        row_spec(seq(1,dim(D)[2]), extra_css = "border-bottom: 2px solid;") %>%
        add_header_above(colnames(nuu), extra_css = border_string) 
      
    }
    
  })
  
  
  #--------------------------CONTINUOUS COVARIATES------------------------------
  output$tablemu <- renderText({
    if (input$choice == "Continuous Covariates Parameters") {
      mu <- data.frame(matrix(ncol = 0,nrow = dim(U)[2]))
      mu$Variable <- paste0("U",1:dim(U)[2])
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
        
        if(C == 3 | C == 4){
          kable_mu <- kable_mu %>% column_spec(4, color = "black",
                               background = my_colors[3])
        }
        if(C == 4){
        kable_mu <- kable_mu %>% column_spec(5, color = "black",
                                             background = my_colors[4])
        }
        
      kable_mu <- kable_mu %>% column_spec(1:(C+1), border_right = "2px solid black", border_left = "2px solid black") %>%
        row_spec(seq(1,dim(U)[2]), extra_css = "border-bottom: 2px solid;") %>%
        column_spec(1, extra_css = "border-left: 2px solid;") %>%
        add_header_above(colnames(mu), extra_css = border_string) 
    }
  })
  
  output$tablesigma <- renderText({
    if (input$choice == "Continuous Covariates Parameters") {
      sigma <- data.frame(matrix(ncol = 0,nrow = dim(U)[2]))
      sigma$Variable <- paste0("U",1:dim(U)[2])
      for(c in 1:C){
        sigma <- cbind(sigma,signif(fitting$params$sigma[,,c],digits = 4),paste0("U",1:dim(U)[2]))
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
      tt <- data.frame(matrix(  rep(c("Variables",paste0("U",1:dim(U)[2])),C) ,ncol = dim(U)[2]*C + C,
                                nrow = 1)  )
      colnames(tt) <- vec
      sigma2 <- rbind(tt,sigma)
      sigma1 <- as.matrix(sigma2)
      colnames(sigma2) <- NULL
      kable_sigma <- sigma2 %>% kbl(align = c("c","c","c","c")) %>%
        kable_styling(font_size = 20) %>% column_spec(seq(1,dim(U)[2]*C + C,dim(U)[2]+1), color = "black",
                                                      background = alpha("grey",0.2), bold = TRUE)  %>%
        column_spec(2:3, color = "black",
        background = my_colors[1]) %>%
        column_spec(5:6, color = "black",
                    background = my_colors[2])
        
        if(C == 3 | C == 4){
          kable_sigma <- kable_sigma %>% column_spec(8:9, color = "black",
                                               background = my_colors[3])
        }
      if(C == 4){
        kable_sigma <- kable_sigma %>% column_spec(11:12, color = "black",
                                             background = my_colors[4])
        }
      
        kable_sigma <- kable_sigma %>% row_spec(1, color = "black",
                 background = alpha("grey",0.2), bold = TRUE) %>%
        column_spec(seq(1,dim(U)[2]*C + C,dim(U)[2]+1), border_right = "2px solid black", border_left = "2px solid black") %>%
        column_spec(dim(U)[2]*C + C, border_right = "2px solid black") %>%
        row_spec(c(1,dim(U)[2]+1), extra_css = "border-bottom: 2px solid;") %>%
        add_header_above(vec2) 
    }
  })
  
  output$plotclust <- renderPlot({
    if (input$choice == "Continuous Covariates Parameters") {
      ggplot(data = data, aes(x = x2, y = x1, color = factor(mclust::map(fitting$z)))) +
        geom_point(shape = 16, size = 2) +
        scale_color_manual(name = 'Cluster',values = c("blue", "green", 
                                                       "deeppink")) +
        labs(x = "U2", y = "U1",title = "Clusters") +
        theme_bw() +
        theme(legend.position = "none",axis.text.x = element_text(size = 12,face = "bold"
                                                                  ),
              plot.title = element_text(size = 16),
              axis.title = element_text(size = 12),
              axis.text.y = element_text(size = 12,face = "bold", hjust = 1.5))
    }
  })
  
  #--------------------------CATEGORICAL COVARIATES------------------------------
  if(v == 1 | v == 2 | v == 3){
  output$tablelam1 <- renderText({
    if (input$choice == "Categorical Covariates Parameters") {
    fir <- data.frame(matrix(ncol = 0,nrow = dim(V1)[2]))
    fir$Category <- 1:dim(V1)[2]
    fir1 <- data.frame(cbind(lambda_11 =  signif(fitting$params$lambda[[1]][,,1],digits = 4),
                            lambda_12 = signif(fitting$params$lambda[[1]][,,2],digits = 4),
                            lambda_13 = signif(fitting$params$lambda[[1]][,,3],digits = 4) ))
    fir <- cbind(fir,fir1)
    colnames(fir) <- c("Category", paste0("\u03BB",pedix[1], pedix[1:C]) )
    fir1 <- as.matrix(fir)
    colnames(fir1) <- NULL
    fir1 %>% kbl(align = c("c","c","c","c")) %>%
      kable_styling(font_size = 20) %>% 
     column_spec(1, color = "black",
                      background = alpha("grey",0.2)) %>%
      column_spec(2, color = "black",
         background = my_colors[1]) %>%
      column_spec(3, color = "black",
                  background = my_colors[2]) %>%
      column_spec(4, color = "black",
                  background = my_colors[3]) %>%
      column_spec(seq(1,C+1), border_right = "2px solid black", border_left = "2px solid black") %>%
      row_spec(seq(1,dim(V1)[2]), extra_css = "border-bottom: 2px solid;") %>%
      add_header_above(colnames(fir), extra_css = border_string) 
    }
  })
  }
  
  if(v == 2 | v == 3){
  output$tablelam2 <- renderText({
    if (input$choice == "Categorical Covariates Parameters") {
      sec <- data.frame(matrix(ncol = 0,nrow = dim(V2)[2]))
      sec$Category <- 1:dim(V2)[2]
      sec1 <- data.frame(cbind(lambda_21 = 
                                signif(fitting$params$lambda[[2]][,,1], digits = 4),
                              lambda_22 = signif(fitting$params$lambda[[2]][,,2], digits = 4),
                     lambda_23 = signif(fitting$params$lambda[[2]][,,3], digits = 4) ))
      sec <- cbind(sec,sec1)
      colnames(sec) <- c("Category", paste0("\u03BB",pedix[2], pedix[1:C]) )
      sec1 <- as.matrix(sec)
      colnames(sec1) <- NULL
      sec1 %>% kbl(align = c("c","c","c","c")) %>%
        kable_styling(font_size = 20) %>%
        column_spec(1, color = "black",
                    background = alpha("grey",0.2)) %>%
        column_spec(2, color = "black",
                                                      background = my_colors[1]) %>%
        column_spec(3, color = "black",
                    background = my_colors[2]) %>%
        column_spec(4, color = "black",
                    background = my_colors[3]) %>%
        column_spec(seq(1,C+1), border_right = "2px solid black", border_left = "2px solid black") %>%
        row_spec(seq(1,dim(V2)[2]), extra_css = "border-bottom: 2px solid;") %>%
        add_header_above(colnames(sec), extra_css = border_string) 
    }
  })
  }
  if(v == 3){
  output$tablelam3 <- renderText({
    if (input$choice == "Categorical Covariates Parameters") {
      thi <- data.frame(matrix(ncol = 0,nrow = dim(V3)[2]))
      thi$Category <- 1:dim(V3)[2]
      thi1 <- data.frame(cbind(lambda_31 = 
                                signif(fitting$params$lambda[[3]][,,1], digits = 4),
                              lambda_32 = signif(fitting$params$lambda[[3]][,,2] , digits = 4),
                              lambda_33 = signif(fitting$params$lambda[[3]][,,3], digits = 4) ))
      thi <- cbind(thi,thi1)
      colnames(thi) <- c("Category", paste0("\u03BB",pedix[3], pedix[1:C]) )
      thi1 <- as.matrix(thi)
      colnames(thi1) <- NULL
      thi1 %>% kbl(align = c("c","c","c","c")) %>%
        kable_styling(font_size = 20) %>% 
        column_spec(1, color = "black",
                    background = alpha("grey",0.2)) %>%
        column_spec(2, color = "black",
                  background = my_colors[1]) %>%
        column_spec(3, color = "black",
                    background = my_colors[2]) %>%
        column_spec(4, color = "black",
                    background = my_colors[3]) %>%
        column_spec(seq(1,C+1), border_right = "2px solid black", border_left = "2px solid black") %>%
        row_spec(seq(1,dim(V3)[2]), extra_css = "border-bottom: 2px solid;") %>%
        add_header_above(colnames(thi), extra_css = border_string) 
    }
  })
  }
  
  output$conditionalUI <- renderUI({
    ui_elements <- list()
    categorical_variables <- paste0("V",1:v)
    for (var in 1:v) {
      ui_elements <- append(ui_elements, list(
        HTML(paste0(
          "<div style='text-align: center; padding: 20px;'>",
          "<div style='border: 2px solid black; background-color: whitesmoke; padding: 10px; display: inline-block; border-radius: 10px;'>",
          "<span style='font-size: 28px; font-weight: bold; font-family: 'Arial', sans-serif; color: #0077B6;'>Categorical Variable V", var, "</span>",
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
        row_spec(seq(1, dim(U)[2] + dim(D)[2] +
                       ifelse(exists("V1"), dim(V1)[2] - 1, 0) +
                       ifelse(exists("V2"), dim(V2)[2] - 1, 0) +
                       ifelse(exists("V3"), dim(V3)[2] - 1, 0)), 
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
        row_spec(seq(1,dim(U)[2] + dim(D)[2] +
                       ifelse(exists("V1"), dim(V1)[2] - 1, 0) +
                       ifelse(exists("V2"), dim(V2)[2] - 1, 0) +
                       ifelse(exists("V3"), dim(V3)[2] - 1, 0)), extra_css = "border-bottom: 2px solid;")  %>%
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
        row_spec(seq(1,dim(U)[2] + dim(D)[2] +
                       ifelse(exists("V1"), dim(V1)[2] - 1, 0) +
                       ifelse(exists("V2"), dim(V2)[2] - 1, 0) +
                       ifelse(exists("V3"), dim(V3)[2] - 1, 0)), extra_css = "border-bottom: 2px solid;") %>%
        add_header_above(c("Variables",colnames(fix1)), extra_css = border_string) 
    }
  })
  
  #--------------------------ISING MODEL----------------------------------------
  output$plot1 <- renderPlot({
    if (input$choice == "Ising Model") {
      lab <- paste0("D",1:dim(D)[2])
      qgraph(fitting$params$int[,,1], fade = FALSE, 
             edge.labels = round(fitting$params$int[,,1], 3), 
             edge.label.cex = 1.5,labels = lab, edge.label.color = "black",
             color = alpha("blue", 0.5),edge.width = 0.7, 
             negCol = "red",posCol = "forestgreen",title = "First cluster Fitted network",
             title.cex = 1.5,edge.label.margin = 0.001,label.font = 2,edge.label.font = 2) 
    }
  })
  output$plot2 <- renderPlot({
    if (input$choice == "Ising Model") {
      lab <- paste0("D",1:dim(D)[2])
      qgraph(fitting$params$int[,,2], fade = FALSE, edge.labels = round(fitting$params$int[,,2], 3), 
             edge.label.cex = 1.5,labels = lab, edge.label.color = "black",color = alpha("green", 0.5),
             edge.width = 0.7, negCol = "red",posCol = "forestgreen",weighted = T,
             title = "Second cluster Fitted network",title.cex =1.5,edge.label.margin = 0.001,
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
             title.cex = 1.5,edge.label.margin = 0.001,label.font = 2,edge.label.font = 2)
    }
  })
  
  #--------------------------RANDOM EFFECTS-------------------------------------
  output$plotre1 <- renderPlot({
    if (input$choice == "Random Effects") {
      re_sim <- REsim(fitting$params$models[[1]])
      re_sim$groupFctr <- " "
      plotREsim(re_sim,labs = TRUE) +
        geom_hline(yintercept = 0, color = I("blue"), size = 1.1) +
        theme_minimal() +
        theme(legend.position = "none",axis.text.x = element_text(size = 12,face = "bold",
                                                                  vjust = 3,angle = 0),
              plot.title = element_text(color = "blue",size = 16, vjust = -2),
              axis.title = element_text(size = 12),
              axis.text.y = element_text(size = 12,face = "bold", hjust = 1.5))
      
    }
  })
  output$plotre2 <- renderPlot({
    if (input$choice == "Random Effects") {
      re_sim <- REsim(fitting$params$models[[2]])
      re_sim$groupFctr <- " "
      plotREsim(re_sim,labs = TRUE) +
        geom_hline(yintercept = 0, color = I("green"), size = 1.1) +
        theme_minimal() +
        theme(legend.position = "none",axis.text.x = element_text(size = 12,face = "bold",
                                                                  vjust = 3,angle = 0),
              plot.title = element_text(color = "green",size = 16, vjust = -2),
              axis.title = element_text(size = 12),
              axis.text.y = element_text(size = 12,face = "bold", hjust = 1.5))
      
    }
  })
  output$plotre3 <- renderPlot({
    if (input$choice == "Random Effects") {
      re_sim <- REsim(fitting$params$models[[3]])
      re_sim$groupFctr <- " "
      plotREsim(re_sim,labs = TRUE) +
        geom_hline(yintercept = 0, color = I("deeppink"), size = 1.1) +
        theme_minimal() +
        theme(legend.position = "none",axis.text.x = element_text(size = 12,face = "bold",
                                                                  vjust = 3,angle = 0),
              plot.title = element_text(color = "deeppink",size = 16, vjust = -2),
              axis.title = element_text(size = 12),
              axis.text.y = element_text(size = 12,face = "bold", hjust = 1.5))
      
    }
  })

}
  

# Run the application
shinyApp(ui = ui, server = server)

