#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)
library(shinythemes)

library(plotly)
library(ggplot2)
library(ggthemes)
library(data.table)
source('https://raw.githubusercontent.com/deepssquared/CUNY.MDS/main/CUNY_DATA_608/Final%20Project/cyto_func.R')


vec.cyto.cols = c('cyto_label',
                  'cell_total',
                  'cp',
                  'ploidy',
                  'clone.count',
  'inv16',
't8.21',
't9.11',
'tv.11',
'add11q',
'del11q',
'add11',
'del11',
'del6',
'del6p',
'del7q',
'del7',
'del7p',
'tri17',
'del17p',
'del5q',
'del5',
'del5p',
'tv.3',
'add3q',
'del3q',
'tv.12',
'del12p',
'del1q',
'add6p',
'add9q',
'add8q',
'del8q',
'add21q',
'add16q',
'tri8',
't6.9',
't6.11',
't3.7',
'inv3',
't1.3',
'del11q',
'BCR.amp',
'del20q',
'add18p',
'del13q',
'add2',
'del3q',
'del9q',
'del11q',
'delY',
't14.18',
't14.16',
'inv12',
'del12p',
'del12q',
'add12q',
't5.10',
'tri1',
'tri1q',
'add1',
'idic21',
'dic5.17',
'tri21',
'tri12',
'MYC.del',
'del8',
'add8q',
'dic7.17',
'mar',
'del16q',
'der1.7',
't12.20',
'add13',
'dic3.17',
'idic20q',
'add5q',
'tv.18',
't9.22',
'abn17p',
'del17',
't15.17',
't8.16')

# Define UI for application that draws a histogram
ui <- navbarPage(
  "karyoParser",
  theme = shinytheme("flatly"),
    # Application title
  
  
  #condition = "input.page != 'plot' & input.page != 'parse'",
  fluidRow(
    h2( img(src = "https://github.com/deepssquared/CUNY.MDS/raw/main/CUNY_DATA_608/Final%20Project/hex.png", width = 80), "Overview"),
    p('karyoParser is an R package designed to parse ISCN-formatted karyotypes in bulk. This work has been previously described in abstract form', a('here.', href = "https://ashpublications.org/blood/article/140/Supplement%201/9119/492067/Quantitative-Cytogenetic-Analysis-with-karyoParser")),
                   h2("Tutorial:"),
                   h4("Parsing [Start Here]:"),
                   p("To start, upload a csv or txt (tab-deliminited) table containing ISCN formatted karyotypes, along with MRNs and Accession Numbers. The column names should be: MRN (Medical Record Number), Accession.No (Identifier specific to the report), Modal_Karyotype (string specific to ISCN). Select the parsing command and the parsing functions will break down the karyotypes into tabular format"), 
                   h4("Plotting:"),
                   p("Next, download the parsed output and upload it for plotting functions"),
                   h4("Sample Data:"),
                   
                   p("Sample datasets are available to test functionality:", a("EHR Sample", href = "https://raw.githubusercontent.com/deepssquared/CUNY.MDS/main/CUNY_DATA_608/Final%20Project/inst/extdata/tutorial_dataset_copy.txt"), "and", a("Parsed Output", href = "https://raw.githubusercontent.com/deepssquared/CUNY.MDS/main/CUNY_DATA_608/Final%20Project/inst/extdata/parsed.csv")),
                   h2("Try it Below!"),
                   
  ),
  
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          
            # sliderInput("bins",
            #             "Number of bins:",
            #             min = 1,
            #             max = 50,
            #             value = 30)
       #   radioButtons("dataset", label = "Choose", choices = list("String", "Upload CSV")),
          fileInput("upload", NULL, buttonLabel = "Upload...", multiple = FALSE, accept = c(".csv", ".txt")),
          radioButtons("page", "Mode", c("EHR Parsing" = "parse", "Plotting" = "plot"), selected = character(0)),
          conditionalPanel(condition = "input.page == 'plot'",
                           selectInput("cyto", "Cytogenetic Plot", vec.cyto.cols)
          ),
        ),

        # Show a plot of the generated distribution
        mainPanel(
          conditionalPanel(condition = "input.page == 'parse'" ,
            DT::dataTableOutput("files"),
        #  downloadButton("download", "Download parsed.csv")
        uiOutput("download_button")),
        
          conditionalPanel(condition = "input.page == 'plot'",
            plotlyOutput("parsed_plot")
          )
       
        
        )

    )
)



# Define server logic required to draw a histogram
server <- function(input, output) {
  source('https://raw.githubusercontent.com/deepssquared/CUNY.MDS/main/CUNY_DATA_608/Final%20Project/cyto_func.R')
  
    # output$distPlot <- renderPlot({
    #     # generate bins based on input$bins from ui.R
    #     x    <- faithful[, 2]
    #     bins <- seq(min(x), max(x), length.out = input$bins + 1)
    # 
    #     # draw the histogram with the specified number of bins
    #     hist(x, breaks = bins, col = 'darkgray', border = 'white',
    #          xlab = 'Waiting time to next eruption (in mins)',
    #          main = 'Histogram of waiting times')
    # })
  data <- reactive({
    req(input$upload)
    
    ext <- tools::file_ext(input$upload$name)
    test = switch(ext,
           csv = read.csv(input$upload$datapath),
           txt = read.delim(input$upload$datapath, sep = '\t'),
           validate("Invalid file; Please upload a .csv or txt file")
    )
    if(input$page == "parse") {
    getCyto.ELN(parseKaryotypewrapper(as.data.table(test), extract.report.text = F, report.text = 'Modal_Karyotype'))
    
    } else {
      as.data.table(test)
    }
  })
  

  output$files <-  DT::renderDataTable({
    req(input$page == "parse")
    datatable(data(),  options = list(scrollX = TRUE))
  })
  

  output$download_button <- renderUI({
    if(!is.null(input$upload) & input$page == "parse" ) {
      downloadButton("download", "Download parsed.csv")
    }
  })
  
  output$download <- downloadHandler(
    filename = function() {
      "parsed.csv"
    },
    content = function(file) {
      write.csv(data(), file, row.names = F)
    }
  )
  
  output$parsed_plot = renderPlotly({

    req(input$page == "plot")
    req(input$upload)
    
    dt = as.data.table(data())
    dt = dt[, .SD, .SDcols = c(input$cyto)]
    setnames(dt, "input.var")
    layout(ggplotly(
      ggplot(data=dt[,.(count = .N), by = "input.var"], aes(x=reorder(input.var, -count), y = count)) + geom_bar(stat = "identity", fill = "navyblue") + labs(x = input$cyto) + theme_minimal()
    )
    )
  })
  
  }

# Run the application 
shinyApp(ui = ui, server = server)
