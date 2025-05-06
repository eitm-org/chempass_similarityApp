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





























