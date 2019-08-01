library(shiny)
library(shinythemes)
library(SNFtool)
library("vegan")
source("functions.R")
ui <- fluidPage(theme = shinytheme("yeti"),
  headerPanel("Integrative Microbiomics",windowTitle="Microbiomics"),
  navbarPage("Integrative Microbiomics",id = "inTabset",
             tabPanel("Biome Submission",value = "panel1",fluidRow(
               column(4,fileInput(inputId = "biome1",label = "Biome: 1"),
                      h4(textOutput(outputId = "biome1_text")),
                      plotOutput(outputId = "biome1_plot"),
                      numericInput("k_biome1", label = h4("Number of species"), value = 10,min=1),
                      actionButton(inputId = "go_biome1",label = "Upload Biome: 1")),
               column(4,fileInput(inputId = "biome2",label = "Biome: 2"),
                      h4(textOutput(outputId = "biome2_text")),
                      plotOutput(outputId = "biome2_plot"),
                      numericInput("k_biome2", label = h4("Number of species"), value = 10,min=1),
                      actionButton(inputId = "go_biome2",label = "Upload Biome: 2")),
               column(4,fileInput(inputId = "biome3",label = "Biome: 3"),
                      h4(textOutput(outputId = "biome3_text")),
                      plotOutput(outputId = "biome3_plot"),
                      numericInput("k_biome3", label = h4("Number of species"), value = 10,min=1),
                      actionButton(inputId = "go_biome3",label = "Upload Biome: 3"))),
               fluidRow(
                 column(2,
                        br(),br(),
                        h5("Click on Next to proceed"),
                        actionButton('jumptot2', 'Next'),offset = 9)
               )
             ),
             tabPanel("Parameter Selection",value = "panel2",fluidRow(
               column(6,
                      h4("Choose a Similarity Measure/Metric"),
                      p("A metric/Similarity measure is a measure or function 
                        that defines distances or similarities between each sample.
                        For example, Bray-Curtis Similarity in 
                        case of ecological or Microbiomic samples ")),
               column(6,
                      h4("Choose a Merging Algorithm/Method"),
                      p("A Merging algorithm is a method through which the datasets
                        are integrated")
                      )),
               fluidRow(
                 hr(),
                 column(6,
                        selectInput("metric","Metric/Similarity Measure :",c("Bray-Curtis")),
                        textOutput("metric_out"),
                        br(),
                        br(),
                        hr(),
                        actionButton(inputId = "merge",label = "Merge")),
                 column(6,
                        selectInput("method","Merging Method :",c("SNF","Weighted SNF")),
                        textOutput("method_out"),
                        hr(),
                        h4("Tuneable Parameters"),
                        p("The default value for these parameters are computed based on your input dataset."),
                        uiOutput('ui'))
               )),
             tabPanel("Cluster Visuvalization",value = "panel3",
                      tags$head(
                        tags$style(
                          HTML(".shiny-notification {
              height: 100px;
              width: 800px;
              position:fixed;
              top: calc(50% - 50px);;
              left: calc(50% - 400px);;
            }"
                               )
                        )
                      ),
                      fluidRow(
                        column(6,
                               h4("Modified Average Silhouette Width"),
                               p("The below table represents the average silhouette width for
                                 different number of clusters.The silhouette value is a measure
                                 of how similar an object is to its own cluster (cohesion)
                                 compared to other clusters (separation). The silhouette 
                                 ranges from −1 to +1, where a high value indicates that the
                                 object is well matched to its own cluster and poorly matched 
                                 to neighboring clusters. Since this is the average sihouette 
                                 width, relying solely on these values for number of cluster
                                 selection is not advisable.")),
                        column(6,
                               h4("Optimal Number of Clusters"),
                               p("Two different methods 1.Best Eigen Gap  2. Rotation cost
                                  are used  to identify the optimal/best number of clusters.
                                ")
                               )
                      ),
                      fluidRow(
                        column(6,tableOutput("k_Table")),
                        column(6,tableOutput("cluster_est1"),
                        column(6,tableOutput("cluster_est2")))
                      ),
                      fluidRow(
                        hr(),
                        column(6,
                               plotOutput(outputId = "merged_biome_plot"),
                               actionButton(inputId = "sim_plot","Similarity Plot"),
                               actionButton(inputId = "log_plot","Log Similarity Plot")
                               ),
                        column(6,
                               numericInput("merged_biome", label = h4("Optimal Number of Clusters"),value=2,min = 2),
                               downloadButton(outputId = "download",label="Downlod Label Files"))
                      )
                      ),
             tabPanel("About",value = "panel4",
                      fluidRow(
                        column(6,h3("Integrative Microbiomics"),
                               p("This Webtool is for integrating microbiomes. 
                                To cite this tool use: link."),
                               br(),
                               p("Please contact the author if you face any problems"),
                               h4("Author"),
                               uiOutput("auth"))
                               ))
             ))
server <- function(input, output, session) {
  in_parm<-reactiveValues(tmp_go_biome3=NULL,tmp_go_biome2=NULL,tmp_go_biome1=NULL)
  print(isolate(in_parm$tmp_go_biome1))
  data1=eventReactive(in_parm$tmp_go_biome1,{
    if (is.null(input$biome1$datapath))
    {return(NULL)}
    else{
      read.csv(input$biome1$datapath,header = TRUE,row.names = 1)
    }
  })
  output$biome1_plot<-renderPlot({bar(data1(),input$k_biome1)})
  output$biome1_text<-renderText({paste("Top",input$k_biome1,"species based on their abundance")})
  data2=eventReactive(in_parm$tmp_go_biome2,{
    if (is.null(input$biome2$datapath))
    {return(NULL)}
    else{
      read.csv(input$biome2$datapath,header = TRUE,row.names = 1)
    }
  })
  output$biome2_plot<-renderPlot({bar(data2(),input$k_biome2)})
  output$biome2_text<-renderText({paste("Top",input$k_biome2,"species based on their abundance")})
  data3=eventReactive(in_parm$tmp_go_biome3,{
    if (is.null(input$biome3$datapath))
    {return(NULL)}
    else{
      read.csv(input$biome3$datapath,header = TRUE,row.names = 1)
    }
  })
  output$biome3_plot<-renderPlot({bar(data3(),input$k_biome3)})
  output$biome3_text<-renderText({paste("Top",input$k_biome3,"species based on their abundance")})
  
output$metric_out<-renderText({
  if(input$metric=="Bray-Curtis"){"Bray-Cutis Similarity
  is a statistic used to quantify the compositional dissimilarity
  between two different sites"}})  
output$method_out<-renderText({
  if(input$method=="SNF"){"Similarity Network Fusion, 
    is a new computational method for data integration. 
    SNF first constructs a sample similarity network for each of 
    the data types and then iteratively integrates these networks using
    a novel network fusion method. Working in the sample network space 
    allows SNF to avoid dealing with different scale, collection bias and 
    noise in different data types. Integrating data in a non-linear fashion 
    allows SNF to take advantage of the common as well as complementary i
    nformation in different data types "}
  else if(input$method=="Weighted SNF"){"Weighted SNF, is a modified version of SNF
    that accepts weightage for each biome, this is necessary as each biome may not be
    sequenced to same depth or may not be of equal quality. This method addresses this 
    issue by weigthing them iteratively at each step of SNF. However, this method can
    only be applied if the number of biomes is 3 or more. In the case of 2 biomes,
    Weighted SNF = SNF"}
    })

output$ui<-renderUI({
if(is.null(input$method))
  return()
  switch(input$method,
       "SNF"=tagList(
        numericInput("K_nn","K Neareast Neighbours",value = round(dim(data1())[1]/10)),
        numericInput("t_iter","Number of Iterations",value=20)),
      "Weighted SNF"=tagList(
        numericInput("K_nn","K Neareast Neighbours",value = round(dim(data1())[1]/10)),
        numericInput("t_iter","Number of Iterations",value=20),
        numericInput("weight1","Weight of the Biome 1",value=length(data1())),
        numericInput("weight2","Weight of the Biome 2",value=length(data2())),
        numericInput("weight3","Weight of the Biome 3",value=length(data3()))
        )
      )
  })

data_merge=eventReactive(input$merge,
                        {
                          validate(need(dim(data1())[1]==dim(data2())[1] & dim(data2())[1]==dim(data3())[1],"No. of patients/samples is not consistent between the biomes. Please Check and rerun!"))
                          validate(need(all(row.names(data1())==row.names(data2()) & row.names(data2()) == row.names(data3())),"The patients/samples are not ordered consistently between the biomes. Please Check and rerun!"))
                          if(input$method=="SNF"){
                          withProgress(message = "Merging Biomes",value = 0,
                                      {merge_snf(list(data1(),data2(),data3()),input$K_nn,input$t_iter)}
                          )}
                          else if(input$method=="Weighted SNF"){
                            withProgress(message = "Merging Biomes",value = 0,
                            {merge_wsnf(list(data1(),data2(),data3()),input$K_nn,input$t_iter,c(input$weight1,input$weight2,input$weight3))}
                            )}
                          })

output$k_Table<-renderTable(
  withProgress(message="Calculating Silhouette",value=0,max_k(data_merge()))
  )

est_tab=reactive(withProgress(message = "Calculating Optimal Clusters",
            as.data.frame(estimateNumberOfClustersGivenGraph(data_merge()),row.names = "Optimal Clusters")
                ))

output$cluster_est1<-renderTable({est_tab()[,1:2]
  })

output$cluster_est2<-renderTable({est_tab()[,3:4]
})

observeEvent(input$sim_plot,{
  output$merged_biome_plot<-renderPlot({
  withProgress(message="Plotting the Biomes",biome_plot(data_merge(),input$merged_biome))
  })
  })

observeEvent(input$log_plot,{
  output$merged_biome_plot<-renderPlot({
    withProgress(message="Plotting the Biomes",biome_log_plot(data_merge(),input$merged_biome))
  })
})

label=reactive({label_create(data_merge(),input$merged_biome)})

output$download<-downloadHandler(
   filename = function() {
     paste('labels-', Sys.Date(), '.csv', sep='')
   },
   content = function(con) {
     write.csv(label(), con)
   },
   contentType = "csv"
 )

#----Moving to different tabs-----
observeEvent(input$jumptot2, {
  in_parm$tmp_go_biome3<- abs(round(rnorm(1)*10))
  in_parm$tmp_go_biome2<- abs(round(rnorm(1)*10))
  in_parm$tmp_go_biome1<- abs(round(rnorm(1)*10))
  updateTabsetPanel(session, "inTabset",selected = "panel2")
})

observeEvent(input$merge, {
  updateTabsetPanel(session, "inTabset",selected = "panel3")
})
#------Linking Action Buttons--------------
observeEvent(input$go_biome1,{in_parm$tmp_go_biome1<- abs(round(rnorm(1)*10))})
observeEvent(input$go_biome2,{in_parm$tmp_go_biome2<- abs(round(rnorm(1)*10))})
observeEvent(input$go_biome3,{in_parm$tmp_go_biome3<- abs(round(rnorm(1)*10))})
#----------------------------------------
url <- a("Jayanth Kumar Narayana (jayanth.kumar@tuta.io)", href="https://jayanthkumarblog.wordpress.com/")
output$auth <- renderUI({tagList(url)})
}
shinyApp(ui = ui, server = server)
#Todo
#Check if the dimensions are the same
#Check if the ordering of the index is same if not reindex