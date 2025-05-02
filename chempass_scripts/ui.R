# ui.R

library(shiny)
library(shinyjs)
library(shinyWidgets)




fluidPage( # fluidpage helps rescale the app display depending on the window its being viewed in
  
  #' adding a logo to the title of the app - might be changed eventually
  tags$head(
    tags$style(HTML("
      .title-logo {
        position: absolute;
        top: 10px;
        right: 20px;
        height: 200px;
      }
      .title-container {
        position: relative;
        padding-top: 10px;
        padding-bottom: 10px;
      }
    "))
  ),
  
  
  div(class = "title-container",
      tags$img(src = "logo.jpg", class = "title-logo")
  ),
  # title of app - will be changed
  titlePanel("CHEMPASS; A Tool for Rapid Chemical Structure Similarity Scoring and Clustering from Common Identifiers."),
  p("Upload your list of compounds in .csv/.txt format either as CID/DTXSIDs/SMILES, select your choice of fingerprint generation, set the clustering threshold, and explore cluster-level statistics and heatmaps."),
  p("Outputs include a .pdf of clusters with similarity scores, molecular properties of compounds, png of heatmap and NMDS plots."),
  helpText("Please refresh app to upload a new file or try another example data"),
  shinyjs::useShinyjs(),
  tags$script(HTML("Shiny.setInputValue('disable_ready', true);")),
  uiOutput("download_ui"),
  # sidepanel layout and buttons
  sidebarLayout(
    sidebarPanel(
      helpText("Use example data to test the app without uploading a file."),
      materialSwitch(inputId = "use_example", label = "Use Example Data", value = FALSE, status = "success"),
      
      conditionalPanel(
        condition = "input.use_example == true",
        p("Processing an arbitrary list of 20 CIDs provided as example data"),
      ),
      
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
        
        
      ),
      
      
      
      
      actionButton("fingerprint_button", "2 - Generate Fingerprints"),
      selectInput("fingerprint_type", "Select fingerprint type:", choices = c("ECFP4", "FCFP4")),
      tags$details(
        tags$summary(HTML('Fingerprint: <span style="color:blue; cursor:pointer;">[?]</span>')),  
        p(strong("What are Morgan fingerprints?")),
        p("Morgan fingerprints encode molecular structures by considering circular atom neighborhoods."),
        tags$ul(
          tags$li(strong("ECFP4:"), " Focuses on atom and bond connectivity (radius 2)."),
          tags$li(strong("FCFP4:"), " Focuses on functional classes of atoms like donor/acceptor roles.")
        ),
        p("Please use the refresh button to change fingerprint type selection and observe the heatmap change accordingly.")
      ),
      
      numericInput("cutoff", "Clustering Cutoff", value = 0.2, step = 0.01),
      tags$details(
        tags$summary(HTML('Clustering Cutoff: <span style="color:blue; cursor:pointer;">[?]</span>')),  
        p(strong("What does Clustering Cutoff mean?")),
        p("This value sets the distance threshold."),
        p("• A ", strong("lower cutoff"), " (e.g., 0.2) forms tighter, smaller clusters."),
        p("• A ", strong("higher cutoff"), " (e.g., 0.6) forms larger, looser clusters."),
        p("Adjust the cutoff depending on how strict you want your clustering to be."),
        p("For chemical structure clustering, a cutoff of 0.2–0.4 is often a good starting point."),
        p("Please use the refresh button to change clustering cutoff and observe the Butina clusters and NMDS plots change accordingly.")
      ),
      actionButton("cluster", "3 - Generate Butina Compound Clusters"),
      
      tags$head(
        tags$style(HTML("a.action-button[disabled] {
        pointer-events: none;
                        opacity: 0;}")))
      
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Heatmap", plotOutput("imageOutput1", width = 1500, height = 1500)),
        tabPanel("Butina clusters", plotOutput("imageOutput2", width = 1000, height = 1000)),
        tabPanel("NMDS Plot", plotOutput("imageOutput3", width = 800, height = 800))
      )
    )
  )
)
