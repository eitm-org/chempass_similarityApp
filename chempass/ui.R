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
  p(strong("Tab 1: Chemical Similarity (single list)")),
  p("Upload your list of compounds in .csv/.txt format either as CID/DTXSIDs/SMILES, select your choice of fingerprint generation, set the clustering threshold, and explore cluster-level statistics and heatmaps."),
  p("zip file for download includes a .pdf of clusters with similarity scores, molecular properties of compounds, high resolution images of heatmap and NMDS plots."),
  p(""),
  p(strong("Tab 2: Compare two lists")),
  p("Upload a reference list of SMILES + a target list of SMILES and for each compound in the target list, get the most structurally similar compound from the reference list."),
  
  
  shinyjs::useShinyjs(),
  tags$script(HTML("Shiny.setInputValue('disable_ready', true);")),
  #uiOutput("download_ui"),
  # sidepanel layout and buttons
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        tabPanel("Chemical Similarity (CS)", 
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
                           tags$td("ID"),
                           tags$td("Name")  # empty column
                         ),
                         tags$tr(
                           tags$td("DTXSID5025659"),
                           tags$td("compound 1")  # empty column
                         ),
                         tags$tr(
                           tags$td("DTXSID5025811"),
                           tags$td("compound 2")
                         ),
                         tags$tr(
                           tags$td("6950"),
                           tags$td("compound 3")
                         ),
                         tags$tr(
                           tags$td("8295"),
                           tags$td("compound 4")
                         ),
                         tags$tr(
                           tags$td("C(C(=O)NCC(=O)O)N"),
                           tags$td("compound 5")
                         ),
                         tags$tr(
                           tags$td("CC(=O)NCCCCN"),
                           tags$td("compound 6")
                         )
                       )
                     ), 
                   ),
                   
                   # File input itself
                   p("          "),
                   fileInput("file_upload", label = NULL),
                   actionButton("process_file", "Check File", style = "color: white; background-color: #A96D4B; border-color: black;"),
                   helpText("Choose either 1a or 1b for further processing"),
                   p("          "),
                   tags$details(
                     tags$summary(HTML('<b>More information for 1a:</b> <span style="color:blue; cursor:pointer;">[?]</span>')),
                     
                     p("-Molecular properties are calculated from SMILES using RDKit"),
                     p("-Input can only be a list of SMILES and/or Names of compounds"),
                     p("-If Names not provided by user; Names follow input order; Compound 1 = SMILES_1, Compound 2 = SMILES_2.")
                     
                   ),
                   p("          "),
                   actionButton("noPubChem", "1a - Faster Clustering (Calculated Properties)", style = "color: white; background-color: #A96D4B; border-color: black;", disabled = T),
                   
                   p("          "),
                   
                   tags$details(
                     tags$summary(HTML('<b>More information for 1b:</b> <span style="color:blue; cursor:pointer;">[?]</span>')),
                     
                     p("This is a slower option; needs to connect with PubChem API"),
                     p("Molecules need to exist in PubChem for this option to work"),
                     p("Molecular properties are downloaded from PubChem")
                     
                   ),
                   p("          "),
                   actionButton("process_file_withPubChem", "1b - Slower Clustering (PubChem Properties)", style = "color: white; background-color: #A96D4B; border-color: black;", disabled = T),
                   p("          "),
                   
                   
                  
                   
                 ),
                 
                 
                 
                 
                 p("          "),
                 tags$details(
                   tags$summary(HTML('<b>Select fingerprint type:</b> <span style="color:blue; cursor:pointer;">[?]</span>')),  
                   p(strong("What are Morgan fingerprints?")),
                   #### add more info on which fingerprint to choose here.
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
        tabPanel("Compare two lists",
                 p(""),
                 materialSwitch(inputId = "example_btn", label = "Load Example Files", value = FALSE, status = "danger"),
                 conditionalPanel(
                   condition = "input.example_btn == true",
                   tags$a(href = "Reference_isomers_enentiomers.csv", "Download reference data",
                          style = "color: #007bff; text-decoration: underline;"),
                   tags$a(href = "Target_isomers_enentiomers.csv", "Download target data",
                          style = "color: #007bff; text-decoration: underline;"),
                 ),
                 
                 conditionalPanel(
                   condition = "input.example_btn == false",
                   fileInput("reference_file_upload", label = "Reference File"),
                   fileInput("target_file_upload", label = "Target File"),
                   
                   p(""),
                   p(strong("Please upload a reference and target list of SMILES.")),
                   p("You can obtain SMILES for your list of compounds from the CS tab --> option 1a --> Download zip"),
                   tags$details(
                     tags$summary(HTML('Upload two lists of SMILES: <span style="color:blue; cursor:pointer;">[?]</span>')),  
                     
                     tags$ul(
                       tags$li(strong("Reference File:"), " List of SMILES"),
                       tags$li(strong("Target File:"), " Compounds (as SMILES) that you want compared to reference list"),
                       tags$li(strong("Output:"), " For each compound in the target list, the output pdf will display the most structurally similar compound from the reference list."),
                     ),
                     p(strong("Reference List:")),
                     tags$table(
                       border = 1,
                       style = "text-align:left; margin-top:20px;",
                       
                       tags$colgroup(
                         tags$col(style = "width:50%"),
                         tags$col(style = "width:50%")
                       ),
                       tags$tbody(
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
                     p("          "),
                     p(strong("Target List:")),
                     tags$table(
                       border = 1,
                       style = "text-align:left; margin-top:20px;",
                       
                       tags$colgroup(
                         tags$col(style = "width:50%"),
                         tags$col(style = "width:50%")
                       ),
                       
                       tags$tbody(
                         tags$tr(
                           tags$td("C(CO)NCCO"),
                           tags$td("")
                         ),
                         tags$tr(
                           tags$td("C1CC(NC1)C(=O)O"),
                           tags$td("")
                         )
                       )
                     ),
                   ),
                   p("          "),
                   
                   
                   
                   helpText("If duplicates are found between reference and target list.\nOnly distinct entries will be kept."),
                   
                 ),
                 
                 
                 actionButton("process_two_files", "1 - Process Files", style = "color: white; background-color: #A96D4B; border-color: black;"),
                 
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
                 selectInput("fingerprint_type_2", "", choices = c("ECFP4", "FCFP4")),
                 actionButton("fingerprint_button_2", "2 - Generate Fingerprints", style = "color: white; background-color: #433e2e; border-color: black;", disabled = T),
                 p("          "),
                 uiOutput("download_ui2")
                 
                 ),
        
        
        
      )
      
      
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Heatmap (CS)", 
                 withSpinner(
                   uiOutput("heatmap_container", width = 1500, height = 1500),
                   type = 6,  
                   color = "#636363",  
                   proxy.height = "150px"  
                   )
                 ),
        tabPanel("Silhouette scores (CS)",
                 withSpinner(
                   uiOutput("silhouette_container", width = 800, height = 800),
                   type = 6,  
                   color = "#636363",  
                   proxy.height = "150px"  
                 )
                 
        ),
        
        tabPanel("Butina Clusters (CS)",
                 withSpinner(
                   uiOutput("cluster_pdf", width = 1500, height = 1500),
                   type = 6,  
                   color = "#636363",  
                   proxy.height = "150px"  
                 )
        ),
        tabPanel("NMDS Plot (CS)",
                 withSpinner(
                   uiOutput("NMDS_container", width = 800, height = 800),
                   type = 6,  
                   color = "#636363",  
                   proxy.height = "150px"  
                 )
        ),
        tabPanel("Two list compare",
                 withSpinner(
                   uiOutput("two_list_pdf", width = 1500, height = 1500),
                   type = 6,  
                   color = "#636363",  
                   proxy.height = "150px"  
                 )
        ),
        
      )
    )
  )
)
