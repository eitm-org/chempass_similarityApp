fluidRow(
  column(
    width = 6,  # adjust width as needed
    
    conditionalPanel(
      condition = "input.use_example == false",
      
      # File input label
      tags$label("Upload file"),
      
      # Instructions placed directly above the browse button
      tags$details(
        tags$summary(HTML('File Upload Instructions: <span style="color:blue; cursor:pointer;">[?]</span>')),  
        p("Please upload a list of compounds in text (.txt) or CSV (.csv) format containing one or more of these identifiers."),
        tags$ul(
          tags$li(
            strong("CID:"), 
            " Compound Identification Number (PubChem). eg: ", code("1983"), " (for Acetaminophen)"
          ),
          tags$li(
            strong("SMILES:"), 
            " Simplified Molecular Input Line Entry System, a line notation for describing the structure of chemical species. eg: ", code("CC(=O)NC1=CC=C(O)C=C1")
          ),
          tags$li(
            strong("DTXSID:"), 
            " DSSTox Substance Identifier used in the Environmental Protection Agency CompTox Dashboard. eg: ", code("DTXSID0020232")
          )
        ),
        
        tags$img(src = "example_table_image.png", width = "70%", height = "auto"), 
      ),
      
      # File input itself
      fileInput("file_upload", label = NULL)
    )
  )
)



conditionalPanel(
  condition = "input.use_example == false",
  fileInput("file_upload", "Upload File", accept = c(".csv", ".txt")),
  actionButton("process_file", "1 - Process This File"),
  
  tags$details(
    tags$summary(HTML('File Upload Instructions: <span style="color:blue; cursor:pointer;">[?]</span>')),  
    p("Please upload a list of compounds in text (.txt) or CSV (.csv) format containing one or more of these identifiers."),
    tags$ul(
      tags$li(
        strong("CID:"), 
        " Compound Identification Number (PubChem). eg: ", code("1983"), " (for Acetaminophen)"
      ),
      tags$li(
        strong("SMILES:"), 
        " Simplified Molecular Input Line Entry System, a line notation for describing the structure of chemical species. eg: ", code("CC(=O)NC1=CC=C(O)C=C1")
      ),
      tags$li(
        strong("DTXSID:"), 
        " DSSTox Substance Identifier used in the Environmental Protection Agency CompTox Dashboard. eg: ", code("DTXSID0020232")
      )
    ),
    
    tags$img(src = "example_table_image.png", width = "70%", height = "auto"), 
  ),
  
  
)

##############################################################
## interactive plots
##############################################################

output$imageOutput4 <- renderPlot({
  req(event3_trigger())
  req(tanimoto_matrix())
  req(clustering_done())
  
  metaMDS_results.points <- metaMDS_results.points_slot()
  
  g <- metaMDS_results.points %>% 
    ggplot(., aes(x = MDS1, y = MDS3)) +
    geom_point(size = 2, col = "grey") +
    geom_point(data = filter(metaMDS_results.points, cluster_ID_toPlot == "plot"), aes(x = MDS1, y = MDS3, col = cluster_ID), size = 4) +
    #scale_color_manual(values = c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#bf5b17','#666666')) +
    theme_classic() +
    theme(text = element_text(size = 15)) +
    labs(title = "NMDS Plot from Distance Matrix",
         x = "MDS1",
         y = "MDS2")
  
  #plot(g)
  
  NMDS_slot(g)
  
})






observeEvent(input$file_upload, {
  
  req(input$file_upload)
  
  file_ext <- tools::file_ext(input$file_upload$name)
  
  
  
  # Try reading the file
  tryCatch({
    mixed_input <- read.csv(input$file_upload$datapath, header = FALSE, stringsAsFactors = FALSE)
    
    # Check number of columns
    if (ncol(mixed_input) != 1) {
      sendSweetAlert(
        session,
        title = "⚠️ Incorrect format",
        text = paste("Expected 1 column, but found ", ncol(mixed_input)),
        type = "warning",
        #timer = 4000,
        #showConfirmButton = FALSE
      )
      return(NULL)
    }
    
    # Clean and name the column
    mixed_input$V1 <- trimws(mixed_input$V1)
    colnames(mixed_input) <- "V1"
    
    # Proceed with your downstream logic here
    # reactive_df(mixed_input) or similar
    
  }, error = function(e) {
    sendSweetAlert(
      session,
      title = "❌ Error reading file",
      text = paste("Please ensure the file is a valid plain-text CSV or TXT.\n\n", e$message),
      type = "error",
      #timer = 6000,
      #showConfirmButton = FALSE
    )
  })
})


library(shinycssloaders)

tabPanel("Heatmap",
         withSpinner(
           plotOutput("imageOutput1", width = 1500, height = 1500),
           type = 6,  # optional: spinner style
           color = "#0dc5c1",  # optional: spinner color
           proxy.height = "150px"  # placeholder space before plot loads
         )
)





library(shiny)

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      h4("Resources"),
      tags$a(
        href = "https://example.com/yourdata.csv", 
        "🔗 Link to XX data", 
        target = "_blank", 
        style = "color: #007bff; text-decoration: underline;"
      )
    ),
    mainPanel(
      # Your main panel content
    )
  )
)

server <- function(input, output, session) {}

shinyApp(ui, server)







ui <- fluidPage(
  checkboxInput("show_table", "Show Table", value = FALSE),
  
  conditionalPanel(
    condition = "input.show_table == true",
    
    tags$table(
      border = 1,
      style = "width:50%; text-align:center; margin-top:20px;",
      
      tags$thead(
        tags$tr(
          tags$th("Column A"),
          tags$th("Column B")
        )
      ),
      tags$tbody(
        tags$tr(
          tags$td("Data 1A"),
          tags$td("Data 1B")
        ),
        tags$tr(
          tags$td("Data 2A"),
          tags$td("Data 2B")
        )
      )
    )
  )
)



if (nrow(dtxsid_df) > 0) {
  pubchem_dtx <- list()
  dtxsids <- dtxsid_df$V1
  
  for (dtxsid in dtxsids) {
    
    cid <- tryCatch({
      get_cid_from_dtxsid(dtxsid)
    }, error = function(e) {
      message(paste("Failed to fetch CID for:", dtxsid, "| Error:", e$message))
      NA  # use NA if failed
    })
    
    pubchem_dtx <- append(pubchem_dtx, list(
      list(input = dtxsid, DTXSID = dtxsid, CID = cid)
    ))
  }
  
  df1 <- pd$DataFrame(pubchem_dtx)
}

cid <- tryCatch({
  get_cid_from_dtxsid(dtxsid)
}, error = function(e) {
  message("Retrying once after error...")
  Sys.sleep(1)
  tryCatch(get_cid_from_dtxsid(dtxsid), error = function(e2) {
    message("Second attempt failed.")
    NA
  })
})









