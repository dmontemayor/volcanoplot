# Rshiny app to visualize differential expression data
# Author: Daniel Montemayor
#


#Requirements
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinybusy)
library(ggplot2)

#random seed
set.seed(44701)

#Set max file size
options(shiny.maxRequestSize=200*1024^2)

#correction methods
corrmethods = c("None", "Bonferroni", "Benjamini-Hochberg")

#init reactive values
#rdata <- reactiveValues(idx = 0, dfagg = NULL, df = NULL, path = NULL)
rdata <- reactiveValues(dfagg = NULL, classes = NULL)

#In-App documentation
doctext<-'
## Summary
1) Visualize differential expression data between 2 classes via volcano plot
2) upload data with samples in rows and with ID in 1st column, class labels in 2nd column
and measured features in the remaining columns. The file should have column names in first row.
3) Select target class and reference class to comapare
4) Enter sigificance threashold
5) Select adjustment for multiple tests: None, Bonferoni, or BH correction
'



#header
header <- dashboardHeader(
    title = "Differential Expression Analysis"
)

#sidebar
sidebar <- dashboardSidebar(
    sidebarMenu(
        menuItem("Upload csv File(s)", tabName = "upload", icon = icon("file-upload")),
        menuItem("Volcano plot", tabName = "analysis", icon = icon("chart-bar")),
        menuItem("Changelog", tabName = "changelog", icon = icon("changelog"))
    )
)

#body
body <-   dashboardBody(
    # Items for each sidebar tab
    tabItems(
        #Process File tab
        tabItem(tabName = "upload",
                fluidRow(
                    #file upload widget
                    column(4,
                           fileInput("file", "Select csv file(s) to process",
                                     multiple = TRUE, accept = ".csv")
                    ),
                    
                    #process data
                    actionBttn("do", "Process file(s)", 
                               style="stretch", icon = icon("calculator"),
                               color="primary", block=TRUE),
                    
                    #message
                    column(12, textOutput("uploadmsg")),
                    
               )
        ),#end upload tabItem
        
        #Analysis tab
        tabItem(tabName = "analysis",
                fluidRow(
                    #Select target class
                    column(4,
                           selectInput(
                               inputId = "targetclass",
                               label = "Select target class of comparison.",
                               choices = NA,
                           )
                    ),
                    #Select reference class
                    column(4,
                           selectInput(
                               inputId = "refclass",
                               label = "Select reference class to compare against.",
                               choices = NA,
                           )
                    ),
                    #Select correction for multiple tests
                    column(4,
                           selectInput(
                               inputId = "corrmethod",
                               label = "Select correction method for multiple tests.",
                               choices = corrmethods,
                           )
                    ),
                    
                    #define significance threashold
                    column(4,
                           numericInput(
                               inputId = "pval",
                                label = "Enter p-value significance cutoff.",
                                value = 0.05,
                            min = 0.,
                            max = 1.,
                        )
                    ),
                    
                    #generate plot
                    actionBttn("plot", "Generate Volcano Plot", 
                               style="stretch", icon = icon("plot"),
                               color="primary", block=TRUE),
                    
                    #download processed data link
                    column(12, downloadButton("downloadData", "Download Processed Data.")),

                    #download processed data link
                    column(12, downloadButton("downloadPlot", "Download Volcano Plot.")),

                    #visualize plot 
                    plotOutput("volcanoplot")
                )
        ),#end analysis tabItem
        
        #Changelog tab
        tabItem(tabName = "changelog", fluidPage( markdown(doctext)))
        
    )#end tabItems
)

#ui
ui <- dashboardPage(header, sidebar, body)


#server
server <- function(input, output, session) {

    #fetch CRPM logo
    output$logo <- renderUI({
        tags$logo(src = "https://dmontemayor.github.io/assets/Long_SOM/horizontal/JPG/UTHSA_Long-SOM_H_CMYK.jpg")
    })
    
    #Process files
    observeEvent(input$do,{
        # show modal spinner to prevent user from interacting while busy
        show_modal_spinner()
        #check if file is loaded
        if(is.null(input$file)){
            #update upload message
            output$uploadmsg <- renderText({"No file Uploaded!"})
        }
        else{
            #init processed files
            outfiles <- c()
            #init classes
            #classes <- c()
            #init features
            #features <- c()
            
            #loop over files
            for (ifile in c(1:length(input$file$datapath))){
                sfile = input$file$datapath[ifile]
                id = input$file$name[ifile]
            
                #get file name
                outfile <- paste(id,".csv", sep="")
                outfiles <- c(outfiles, outfile)
                
                #get file data
                filepath <- sfile
                data <- read.csv(file = filepath)
                
                #assert features values are numeric
                fslice <- 3:ncol(data)
                data[fslice] <- lapply(data[fslice], as.numeric)
                
                #assert class column (column 2) is named "class"
                names(data)[2] <- "class"

                #accumulate output data    
                if(is.null(rdata$dfagg)){
                    rdata$dfagg <- data
                }
                else{
                    rdata$dfagg <- merge(rdata$dfagg, data, by="class")
                }
            }
            
            #get accumulated classes
            classes = unique(rdata$dfagg["class"])
            
            #get accumulated features
            rdata$features = unique(names(rdata$dfagg)[3:ncol(rdata$dfagg)])
            
            #update target classes for ui 
            updateSelectInput(session, "targetclass", "Select target class of comparison.", choices = classes)
            #update ref classes for ui 
            updateSelectInput(session, "refclass", "Select reference class to compare against.", choices = classes)
            
            
            #update upload message
            output$uploadmsg <- renderText({paste(c(length(outfiles), "Files processed:", paste(outfiles, collapse = ", ")), collapse = " ")}) 

            ## Downloadable csv of aggregated dataset ----
            #output$downloadData <- downloadHandler(
            #    filename = function() {
            #        file.path(getwd(), "volcanoplot_data.csv")
            #    },
            #    content = function(file) {
            #        write.csv(rdata$dfagg, file, row.names = FALSE)
            #    }
            #)
           
        }
        #remove spinner when done
        remove_modal_spinner()
    })#end do proccess file
    
            
    #Volcano plot
    observeEvent(input$plot,{

        #init dataframe
        results <- data.frame(matrix(NA, nrow=length(rdata$features), ncol=8))
        names(results) <- c("target_mean","ref_mean", "difference", "CI.lower", "CI.upper", "p.value", "log2FoldChange", "neglog10p")
        row.names(results) <- rdata$features
        
        #loop over features
        for(feat in rdata$features){
            x <- rdata$dfagg[rdata$dfagg$class==input$targetclass, feat]
            y <- rdata$dfagg[rdata$dfagg$class==input$refclass, feat]
            z<- t.test(x, y, alternative = "two.sided", var.equal = FALSE)
            results[feat, "target_mean"] <- z$estimate[1] 
            results[feat, "ref_mean"] <- z$estimate[2] 
            results[feat, "difference"] <- z$estimate[1] - z$estimate[2]
            results[feat, "CI.lower"] <- z$conf.int[1]
            results[feat, "CI.upper"] <- z$conf.int[2]
            results[feat, "p.value"] <- z$p.value
            results[feat, "log2FoldChange"] <- log2(z$estimate[1] / z$estimate[2])
            results[feat, "neglog10p"] <- -log10(z$p.value)
        }
        
        #adjust for multiple test
        #if(input$cormethod==corrmethods[1]){
        #    #do nothing
        #}
        #if(input$cormethod==corrmethods[2]){
        #    #Bonferoni
        #}
        #if(input$cormethod==corrmethods[3]){
        #    #BH
        #}
        
        # Downloadable processed data ----
        output$downloadData <- downloadHandler(
            filename = function() { paste("volcanoplot", '.csv', sep='') },
            content = function(file) {
                write.csv(results, file)
            }
        )
    
        #generate plot
        plotInput <- reactive({
            #volcano plot
            volcano = ggplot(data = results, aes(x = log2FoldChange, y = neglog10p))
            volcano + geom_point() + 
                geom_text(aes(label=row.names(results)),hjust=0, vjust=0) + 
                geom_hline(yintercept=-log10(input$pval), linetype="dashed")+ 
                geom_vline(xintercept=1, linetype="dashed")+
                geom_vline(xintercept=-1, linetype="dashed")+
                ggtitle("Atrophic vs Normal Tubules") +
                xlab("Log2 Fold Change") + ylab("- log10 pvalue")
        })
        
        output$volcanoplot <- renderPlot({
            print(plotInput())
        })
        
        # Downloadable volcano plot ----
        output$downloadPlot <- downloadHandler(
            filename = function() { paste("volcanoplot", '.png', sep='') },
            content = function(file) {
                device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
                ggsave(file, plot = plotInput(), device = device)
            }
        )
                
    })#end plot 
 
                

}#end Server

# Run the application 
shinyApp(ui = ui, server = server)
