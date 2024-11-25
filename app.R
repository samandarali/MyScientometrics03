# List of required libraries

## Aria, M. & Cuccurullo, C. (2017) bibliometrix: An R-tool for comprehensive science mapping analysis, Journal of Informetrics, 11(4), pp 959-975, Elsevier.


if(!require(shiny)){
  install.packages("shiny")
  library(shiny)
}

if(!require(bibliometrix)){
  install.packages("bibliometrix")
  library(bibliometrix)
}

if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}

if(!require(RColorBrewer)){
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

if(!require(ggthemes)){
  install.packages("ggthemes")
  library(ggthemes)
}

if(!require(gridExtra)){
  install.packages("gridExtra")
  library(gridExtra)
}

if(!require(grid)){
  install.packages("grid")
  library(grid)
}

if(!require(reshape2)){
  install.packages("reshape2")
  library(reshape2)
}

if(!require(treemapify)){
  install.packages("treemapify")
  library(treemapify)
}

if(!require(shinylogs)){
  install.packages("shinylogs")
  library(shinylogs)
}

ui <- fluidPage(
  use_tracking(),
  titlePanel("MyScientometrics"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose a file", accept = c(".txt")),
      br(),
      helpText("Please upload a Web of Science file in plain text format."),
      br(),
      actionButton("biblioButton", "Run Bibliometric Analysis"),
      hr(),
      helpText("The bibliometric analysis may take several minutes to run. Please do not close or refresh the page until the analysis is complete."),
      hr(),
      helpText("The analysis results will be displayed below."),
      checkboxInput("use_sample", "Use Sample Dataset", value = FALSE)
      
      
    ),
    
    mainPanel(
      plotOutput("trendPlot", height = "300px", width = "500px"),
      plotOutput("WoS", height = "200px", width = "700px"),
      plotOutput("country", height = "200px", width = "500px"),
      plotOutput("journal", height = "200px", width = "500px"),
      plotOutput("Pro.Authors", height = "200px", width = "300px"),
      plotOutput("Net.Authors", height = "300px", width = "500px"),
      plotOutput("Net.keywords", height = "300px", width = "500px")
      
      
    )
  )
)

server <- function(input, output) {
  options(shiny.maxRequestSize = 50 * 1024^2)  
  track_usage(storage_mode = store_json(path = "logs/"))
  
  # Load a sample dataset
  sample_data_path <- "PatentAnalysis.txt"  
  
  
  
  
  # Convert data to bibliometrix format
  biblio_data <- reactive({
    # Use sample dataset if no file is uploaded
    file_path <- if (!is.null(input$file) && !input$use_sample) {
      input$file$datapath
    } else {
      sample_data_path
    }
    
    
    df <- convert2df(file_path, dbsource = "wos", format = "plaintext")
    
    
    
    #    df <- convert2df(input$file$datapath, dbsource = "wos", format = "plaintext")
    df <- subset(df, PY < 2023)
    results <- biblioAnalysis(df, sep = ";")
    summary_results <- summary(results, k = 10, pause = FALSE)
    # preparing data for trend graph
    trend_data <- data.frame(Year = as.integer(gsub("\\s+", "", summary_results$AnnualProduction$Year)),
                             Articles = summary_results$AnnualProduction$Articles)
    TcitationY <- df[,c("PY", "TC")]
    Citation <- aggregate(TcitationY$TC, by=list(TcitationY$PY), FUN=sum)
    colnames(Citation) <- c("Year", "Citation")
    aggregate1 <- merge(Citation, trend_data, by.Citation="Year")
    aggregate1$T.Citation <- (aggregate1$Citation * aggregate1$Articles)
    
    #############
    country.citation<-data.frame("country"=summary_results$MostProdCountries$Country,
                                 "Articles"=as.numeric(summary_results$MostProdCountries$Articles),
                                 "SCP" =as.numeric(summary_results$MostProdCountries$SCP),
                                 "MCP" =as.numeric(summary_results$MostProdCountries$MCP)
    )
    
    df_long <- tidyr::gather(country.citation, key = "Variable", value = "Value", -country, -Articles)
    df_long$country<-as.factor(df_long$country)
    df_long$Variable<-as.factor(df_long$Variable)
    ###############
    
    
    # preparing data for Treemap
    
    df.WOS <- df$SO
    
    # split the `WoS Categories` column into `Topic`
    df.WOS <- data.frame(df.WOS)
    categories <- unlist(strsplit(as.character(df.WOS$df.WOS), "; "))
    # count the occurrence of each category
    category_count <- data.frame(table(categories))
    
    ######
    category_count2<-category_count[order(category_count$Freq, decreasing = TRUE), ] #ordering by frequency
    n.category=6
    Wos.cat.short<-head(category_count2,  n.category=6)
    
    others.count<-nrow(category_count)-n.category
    
    # Calculate the sum of Freq for all rows
    total_freq <- sum(category_count2$Freq)
    
    # Calculate the sum of Freq for the top 5 rows
    top5_freq <- sum(Wos.cat.short$Freq)
    
    # Calculate the sum of Freq for the rest of the rows
    rest_freq.av <- (total_freq - top5_freq)/others.count
    # Create a new row for the rest_freq value
    rest_row <- data.frame(categories = "average Other fields", Freq = rest_freq.av)
    
    # Add the rest_row to the category_count2 dataframe
    Wos.cat.short <- rbind(Wos.cat.short, rest_row)
    
    ######Top journal of the field
    journal<-summary_results$MostRelSources
    journal$`Sources       `<-as.factor(journal$`Sources       `)
    journal$Articles<-as.numeric(journal$Articles)
    
    ###### Most Productive Authors
    Authors.pro<-summary_results$MostProdAuthors
    colnames(Authors.pro)<-c("Authors", "Articles", "S.Authors", "Articles Fractionalized")
    Authors.pro$Authors<-as.factor(Authors.pro$Authors)
    
    Authors.pro$Articles<-as.numeric(Authors.pro$Articles)
    
    ##### AUTHORS COLLABRATION NETWORK
    NetMatrix.authors <- biblioNetwork(df, analysis = "collaboration", 
                                       network = "authors", sep = ";")
    
    ################# Keywords Co-occurrences
    
    NetMatrix.keywords <- biblioNetwork(df, analysis = "co-occurrences", network = "keywords", sep = ";")
    
    
    ## return funcgtion for server
    return(list(aggregate1 = aggregate1, df_long=df_long, Wos.cat.short=Wos.cat.short, journal=journal, Authors.pro=Authors.pro , NetMatrix.authors=NetMatrix.authors , NetMatrix.keywords=NetMatrix.keywords
    ) )
    
  })
  
  # Plot the trend graph
  output$trendPlot <- renderPlot({
    if(is.null(biblio_data()$aggregate1)) {
      return(NULL)
    }
    coeff <- 50
    ggplot(biblio_data()$aggregate1, aes(x = Year)) +
      geom_area(aes(y = Articles), fill = "#69b3a2", alpha = 0.5) +
      geom_point(aes(y = Articles), color = "#69b3a2") +
      geom_line(aes(y = Articles), color = "#69b3a2") + 
      geom_line(aes(y = Citation / coeff), color = "deeppink3") +
      scale_y_continuous(
        name = "Number of Publication",
        sec.axis = sec_axis(~.*coeff, name = "Number of Citation")
      ) +
      ggtitle("Number of publication and their citation over time") +
      theme(legend.text = element_text(size = 10),
            axis.title.y = element_text(color = "#69b3a2", size = 13),
            axis.title.y.right = element_text(color = "deeppink3", size = 13),
            axis.text.x = element_text(angle = 45, hjust = 1), 
            panel.background = element_blank(),
            panel.grid.major = element_line(colour = "gray95", size = 0.2, linetype = "solid"),
            panel.border = element_rect(color = "gray", fill = NA)) +
      scale_x_continuous(breaks = seq(2010, 2022, by = 1))
  })
  
  
  # Plot the WoS categories
  output$WoS <- renderPlot({
    if(is.null(biblio_data()$Wos.cat.short)) {
      return(NULL)
    }
    # Treemap
    
    colourCount = 11
    getPalette = colorRampPalette(brewer.pal(9, "Set3"))
    
    Tree_chart<-ggplot(biblio_data()$Wos.cat.short, aes(area = Freq, fill = categories, label = paste(categories, Freq, sep = "\n"))) +
      geom_treemap() +
      geom_treemap_text(colour = "black",
                        place = "centre",
                        size = 15)+ 
      scale_fill_manual(values = getPalette(colourCount))+
      # scale_fill_brewer(palette = "BrBG")+
      ggtitle("articles by Wos Categories")+
      theme(legend.position = "none", legend.text=element_text(size=10))
    
    
    
    Tree_chart
  })
  
  ### Number of articles by country
  
  output$country <- renderPlot({
    if(is.null(biblio_data()$df_long)) {
      return(NULL)
    }
    coeff <- 50
    ggplot(biblio_data()$df_long, aes(x =reorder(country, -Value) , y =Value , fill = Variable),) +
      geom_bar( stat = "identity", position = "stack") +
      scale_fill_manual(values = c("hotpink3", "cyan4"), 
                        name = "Variable",
                        labels = c("SCP", "MCP")) +
      labs(x = "country", y = "Number of Articles", title = "Most Productive Countries") +
      theme_minimal()+
      coord_flip()+
      geom_text(x = Inf, y = Inf, label = "SCP: Single Country Publications, MC: Multiple Country Publications", 
                hjust = 1, vjust = 1, size = 3, color="darkgoldenrod4")
    
    
    
  })
  
  ######Top journal of the field
  
  output$journal <- renderPlot({
    if(is.null(biblio_data()$journal)) {
      return(NULL)
    }
    ggplot(data=biblio_data()$journal, aes(y=Articles, x=reorder(`Sources       `, -Articles))) +
      geom_bar(stat="identity", fill="cyan4")+
      labs(title="Top journal of the field", x="", y="")+
      coord_flip()+
      theme_minimal()
    
  })
  
  ###### Most Productive Authors
  output$Pro.Authors <- renderPlot({
    if(is.null(biblio_data()$Authors.pro)) {
      return(NULL)
    }
    ggplot(data=biblio_data()$Authors.pro, aes(y=Articles, x=reorder(Authors, -Articles))) +
      geom_bar(stat="identity", fill="cyan4")+
      labs(title="Most Productive Authors", x="", y="")+
      coord_flip()+
      theme_minimal()
    
  })
  
  # Authors collarobation
  output$Net.Authors <- renderPlot({
    if(is.null(biblio_data()$NetMatrix.authors)) {
      return(NULL)
    }
    # Treemap
    
    networkPlot(biblio_data()$NetMatrix.authors, n = 30, type = "kamada", Title = "Authors Collaboration",
                #                    label = FALSE,
                halo=FALSE,
                labelsize=0.7
    ) 
    
  })
  
  output$Net.keywords <- renderPlot({
    if(is.null(biblio_data()$NetMatrix.keywords)) {
      return(NULL)
    }
    # Treemap
    
    networkPlot(biblio_data()$NetMatrix.keywords, normalize="association", weighted=T, n = 30, Title = "Keyword Co-occurrences", 
                type = "fruchterman", size=T,edgesize = 5,labelsize=0.7)
  })  
  
  
}

shinyApp(ui, server)