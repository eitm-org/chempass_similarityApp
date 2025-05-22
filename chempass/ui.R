# ui.R

library(shiny)
library(shinyjs)
library(shinyWidgets)
library(plotly)
library(shinybusy)
library(shinycssloaders)


fluidPage( # fluidpage helps rescale the app display depending on the window its being viewed in
  
  #' adding a logo to the title of the app - might be changed eventually
  tags$head(
    tags$style(HTML("
      .title-logo {
        position: absolute;
        top: 10px;
        right: 20px;
        height: 150px;
      }
      .title-container {
        position: relative;
        padding-top: 5px;
        padding-bottom: 5px;
      }
    "))
  ),
  
  
  div(class = "title-container",
      tags$img(src = "logo.jpg", class = "title-logo")
  ),
  # title of app - will be changed
  titlePanel("CHEMPASS; A Tool for Chemical Structure Similarity Scoring and Clustering."),
  p("Upload your list of compounds in .csv/.txt format either as CID/DTXSIDs/SMILES, select your choice of fingerprint generation, set the clustering threshold, and explore cluster-level statistics and heatmaps."),
  p("zip file for download includes a .pdf of clusters with similarity scores, molecular properties of compounds, high resolution images of heatmap and NMDS plots."),
  helpText("Input data may be stored temporarily for testing and quality improvement purposes only."),
  shinyjs::useShinyjs(),
  tags$script(HTML("Shiny.setInputValue('disable_ready', true);")),
  #uiOutput("download_ui"),
  # sidepanel layout and buttons
  sidebarLayout(
    sidebarPanel(
      helpText("Use example data to test the app without uploading a file."),
      materialSwitch(inputId = "use_example", label = "Use Example Data", value = FALSE, status = "danger"),
      
      conditionalPanel(
        condition = "input.use_example == true",
        p("Processing an arbitrary list of 15 CIDs provided as example data"),
        tags$a(href = "CID_example_data.csv", "Download example data",
               style = "color: #007bff; text-decoration: underline;"),
      ),
      
      conditionalPanel(
        condition = "input.use_example == false",
        
        # File input label
        #tags$label("Upload File"),
        
        # Instructions placed directly above the browse button
        tags$details(
          tags$summary(HTML('<b>Upload File:</b> <span style="color:blue; cursor:pointer;">[?]</span>')),  
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
          
          tags$table(
            border = 1,
            style = "text-align:left; margin-top:20px;",
            
            tags$colgroup(
              tags$col(style = "width:50%"),
              tags$col(style = "width:50%")
            ),
            
            tags$tbody(
              tags$tr(
                tags$td("DTXSID5025659"),
                tags$td("")  # empty column
              ),
              tags$tr(
                tags$td("DTXSID5025811"),
                tags$td("")
              ),
              tags$tr(
                tags$td("6950"),
                tags$td("")
              ),
              tags$tr(
                tags$td("8295"),
                tags$td("")
              ),
              tags$tr(
                tags$td("C(C(=O)NCC(=O)O)N"),
                tags$td("")
              ),
              tags$tr(
                tags$td("CC(=O)NCCCCN"),
                tags$td("")
              )
            )
          ), 
        ),
        
        # File input itself
        fileInput("file_upload", label = NULL),
        actionButton("process_file", "1 - Process This File", style = "color: white; background-color: #A96D4B; border-color: black;")
      ),
      
     
      
      
      p("          "),
      tags$details(
        tags$summary(HTML('<b>Select fingerprint type:</b> <span style="color:blue; cursor:pointer;">[?]</span>')),  
        p(strong("What are Morgan fingerprints?")),
        p("Morgan fingerprints encode molecular structures by considering circular atom neighborhoods."),
        tags$ul(
          tags$li(strong("ECFP4:"), " Focuses on atom and bond connectivity (radius 2)."),
          tags$li(strong("FCFP4:"), " Focuses on functional classes of atoms like donor/acceptor roles.")
        ),
      ),
      selectInput("fingerprint_type", "", choices = c("ECFP4", "FCFP4")),
      actionButton("fingerprint_button", "2 - Generate Fingerprints", style = "color: white; background-color: #433e2e; border-color: black;"),
      
      p("          "),
      tags$details(
        tags$summary(HTML('<b>Clustering Cutoff:</b> <span style="color:blue; cursor:pointer;">[?]</span>')),  
        p(strong("What does Clustering Cutoff mean?")),
        p("This value sets the distance threshold."),
        p("• A ", strong("lower cutoff"), " (e.g., 0.2) forms tighter, smaller clusters."),
        p("• A ", strong("higher cutoff"), " (e.g., 0.6) forms larger, looser clusters."),
        p("Adjust the cutoff depending on how strict you want your clustering to be."),
        p("For chemical structure clustering, a cutoff of 0.2–0.4 is often a good starting point."),
        p("For smaller datasets with fewer than 20 compounds, you may need to increase the cutoff beyond 0.4")
      ),
      numericInput("cutoff", "", value = 0.5, step = 0.01),
      actionButton("cluster", "3 - Generate Butina Clusters", style = "color: white; background-color: #3a7e9b; border-color: black;"),
      p("          "),
      uiOutput("download_ui"),
      
      tags$head(
        tags$style(HTML("a.action-button[disabled] {
        pointer-events: none;
                        opacity: 0;}")))
      
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Heatmap", 
                 withSpinner(
                   uiOutput("heatmap_container", width = 1500, height = 1500),
                   type = 6,  
                   color = "#636363",  
                   proxy.height = "150px"  
                   )
                 ),
        
        tabPanel("Butina Clusters",
                 withSpinner(
                   uiOutput("cluster_pdf", width = 1500, height = 1500),
                   type = 6,  
                   color = "#636363",  
                   proxy.height = "150px"  
                 )
        ),
        tabPanel("NMDS Plot",
                 withSpinner(
                   uiOutput("NMDS_container", width = 800, height = 800),
                   type = 6,  
                   color = "#636363",  
                   proxy.height = "150px"  
                 )
        ),
        
      )
    )
  )
)
