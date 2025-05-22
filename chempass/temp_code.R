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


###############################################################################
### testing chunks of code
###############################################################################

mixed_input <- read.csv("/Users/kdabke/Documents/GitHub/CHEMPASS_Rshiny/test_files/SMILES_test_file3.csv",
                        header = FALSE, stringsAsFactors = FALSE)

mixed_input$V1 <- trimws(mixed_input$V1)


duplicate_rows <- mixed_input[duplicated(mixed_input),]
mixed_input <- mixed_input %>% count(V1) %>% filter(n == 1) %>% subset(select = V1)



numbers_only <- function(x) !grepl("\\D", x)

mixed_input$DTXSID <- grepl("DTXSID", mixed_input$V1)

mixed_input$CID <- numbers_only(mixed_input$V1)

dtxsid_df <- filter(mixed_input, DTXSID == TRUE)
dtxsid_df <- subset(dtxsid_df, select = V1)


cid_df <- filter(mixed_input, CID == TRUE)
cid_df <- subset(cid_df, select = V1)


smiles_df <- filter(mixed_input, DTXSID == FALSE & CID == FALSE)
smiles_df <- subset(smiles_df, select = V1)


if (dim(dtxsid_df)[1] > 0) {
  pubchem_dtx <- list()
  dtxsids <- dtxsid_df$V1
  
  for (dtxsid in dtxsids) {
    cid <- tryCatch({
      get_cid_from_dtxsid(dtxsid)
    }, error = function(e) {
      #message("Retrying once after error...")
      Sys.sleep(1)
      tryCatch(get_cid_from_dtxsid(dtxsid), error = function(e2) {
        #message("Second attempt failed.")
        "error"
      })
    })
    
    pubchem_dtx <- append(pubchem_dtx, list(list(input = dtxsid,'DTXSID' = dtxsid, 
                                                 'CID' = cid, 'TITLE' = NULL)))
    
  }
  
  df1 <- pd$DataFrame(pubchem_dtx)
  df1 <- df1 %>% count(CID) %>% filter(n == 1)
}




if (dim(smiles_df)[1] > 0) {
  pubchem_smiles <- list()
  smiles <- smiles_df$V1
  
  for (smile in smiles) {
    encoded <- base64encode(charToRaw(smile))
    #cid <- get_CID_from_SMILES(encoded)
    
    cid <- tryCatch({
      get_CID_from_SMILES(encoded)
    }, error = function(e) {
      message("Retrying once after error...")
      Sys.sleep(1)
      tryCatch(get_CID_from_SMILES(encoded), error = function(e2) {
        message("Second attempt failed.")
        "error"
      })
    })
    
    pubchem_smiles <- append(pubchem_smiles, list(list(input = smile,'DTXSID' = NULL, 
                                                       'CID' = cid, 'TITLE' = NULL)))
    
  }
  
  df2 <- pd$DataFrame(pubchem_smiles)
  df2 <- df2 %>% count(CID) %>% filter(n == 1)
}



if (dim(cid_df)[1] > 0) {
  pubchem_cid <- list()
  cids <- cid_df$V1
  
  for (cid in cids) {
    title <- check_cid(cid)
    pubchem_cid <- append(pubchem_cid, list(list(input = cid, 'DTXSID' = NULL, 'CID' = cid,
                                                 'TITLE' = title)))
    
  }
  
  df3 <- pd$DataFrame(pubchem_cid)
  df3 <- df3 %>% count(CID) %>% filter(n == 1)
}





dfs_list <- list()
if (exists("df1")) dfs_list$df1 <- df1
if (exists("df2")) dfs_list$df2 <- df2
if (exists("df3")) dfs_list$df3 <- df3

df <- bind_rows(dfs_list)
error_rows <- df %>%
  filter(if_any(everything(), ~ . == 'error'))


df_filt <- filter(df, !(input %in% error_rows$input))
df_filt <- subset(df_filt, select = -c(TITLE))
df_filt <- mutate(df_filt, CID = as.character(CID))

property_df <- get_properties_from_CIDs(paste(df_filt$CID, collapse = ","))
property_df <- mutate(property_df, CID = as.character(CID))


missing_names <- (property_df[which(is.na(property_df$Title)),])$CID


if(length(missing_names) > 0) {
  property_df <- property_df[-c(which(is.na(property_df$Title))),]
  
  
}
synonym_list <- list()
for (cid in df_filt$CID) {
  synonyms <- get_synonyms_from_CID(cid)
  
  synonym_list <- append(synonym_list, list(list('CID' = cid,'SYNONYMS' = synonyms)))
  
}
syn_df <- pd$DataFrame(synonym_list)
syn_df <- mutate(syn_df, CID = as.character(CID))

#print(df_filt)
df_filt <- left_join(df_filt, property_df)
df_filt <- left_join(df_filt, syn_df)


if(length(missing_names) > 0) {
  df_filt <- df_filt[-c(which(is.na(df_filt$Title))),]
}



#print(df_filt)

merged_df <- df_filt
#print(merged_df)
colnames(merged_df)[grep("Title", colnames(merged_df))] <- "Name"
#print(merged_df)

missing_data(c(error_rows$input, missing_names))

merged_df$Name <- gsub(" ", "_", merged_df$Name)
merged_df$ROMol <- lapply(merged_df$SMILES, function(smiles) {
  
  return(Chem$MolFromSmiles(smiles))
})

merged_df <- as.data.frame(merged_df)


###############################################################################
### SMILES not working
###############################################################################

smile <- gsub("\\s+", "", "CN1CC(=O)N=C1N")
encoded <- base64encode(charToRaw(smile))
get_CID_from_SMILES(encoded)


mixed_input <- read.csv("/Users/kdabke/Downloads/CHEMPASS-test_og.csv", header = FALSE, stringsAsFactors = FALSE)

mixed_input$DTXSID <- grepl("DTXSID", mixed_input$V1)

numbers_only <- function(x) !grepl("\\D", x)
mixed_input$CID <- numbers_only(mixed_input$V1)

dtxsid_df <- filter(mixed_input, DTXSID == TRUE)
dtxsid_df <- subset(dtxsid_df, select = V1)


cid_df <- filter(mixed_input, CID == TRUE)
cid_df <- subset(cid_df, select = V1)


smiles_df <- filter(mixed_input, DTXSID == FALSE & CID == FALSE)
smiles_df <- subset(smiles_df, select = V1)


if (dim(dtxsid_df)[1] > 0){
  df1 <- dtxsid_pubchem(dtxsid_df)
}



#' this could be a smiles function
#' usage eg:
#' df2 <- smiles_pubchem(smiles_df)

if (dim(smiles_df)[1] > 0){
  df2 <- smiles_pubchem(smiles_df)
  
}



#' this would be a CID function
#' usage example
#' df3 <- input_cid_reformating(cid_df)

if (dim(cid_df)[1] > 0){
  df3 <- input_cid_reformating(cid_df)
  
}


dfs_list <- list()
if (exists("df1")) dfs_list$df1 <- df1
if (exists("df2")) dfs_list$df2 <- df2
if (exists("df3")) dfs_list$df3 <- df3

#' this could also be a function in itself
#' usage example
#' merged_df <- cid_processing(dfs_list)

df <- bind_rows(dfs_list)

error_rows <- df %>%
  filter(if_any(everything(), ~ . == 'error'))

df_filt <- filter(df, !(input %in% error_rows$input))

df_filt <- subset(df_filt, select = -c(TITLE))


#' mutates the cid into a character

df_filt <- mutate(df_filt, CID = as.character(CID))

df_filt <- df_filt %>% distinct(CID, .keep_all = TRUE)

# get properties from Pubchem
property_df <- tryCatch({
  get_properties_from_CIDs(paste(df_filt$CID, collapse = ","))
}, error = function(e) {
  #message("Retrying once after error...")
  Sys.sleep(1)
  tryCatch(get_properties_from_CIDs(paste(df_filt$CID, collapse = ",")), error = function(e2) {
    #message("Second attempt failed.")
    "error"
  })
})

#' re-enforces that the CIDs is a character column
#' this is to ensure left join works downstream

property_df <- mutate(property_df, CID = as.character(CID))


syn_df <- cid_to_synonym(df_filt)

df_filt <- left_join(df_filt, property_df) # only common column here should be CID
df_filt <- left_join(df_filt, syn_df) # only common column here is CID


#' now to catch the missing names from titles
#' so this checks whether a CID has a title missing
#' and collects those missing into a new variable called missing_names
missing_names <- (df_filt[which(is.na(df_filt$Title)),])$input

#' if there are missing_names, it removes the entire row from the datafrmae
#' this ensures that downstream, the CIDs that get processed are clean

if(length(missing_names) > 0) {
  df_filt <- df_filt[-c(which(is.na(df_filt$Title))),]
}

#' checks the data frame for duplicate CIDs after all the processing


df_filt <- df_filt %>% distinct(CID, .keep_all = TRUE)

merged_df <- df_filt

colnames(merged_df)[grep("Title", colnames(merged_df))] <- "Name"


#' saves the error and missing CIDs and names
#' is it possible to show here the missing input instead
#' might be better for user


#' is this needed if I use clean names somwhere else
#' check the difference here
merged_df$Name <- gsub(" ", "_", merged_df$Name)
merged_df$ROMol <- lapply(merged_df$SMILES, function(smiles) {
  
  return(Chem$MolFromSmiles(smiles))
})

merged_df <- as.data.frame(merged_df)

fpgen <- AllChem$GetMorganGenerator(radius = as.integer(2), fpSize = as.integer(2048))
fingerprints <- list()
for (i in 1:nrow(merged_df)){ 
  mol <- merged_df$ROMol[i]
  fingerprints[[i]] <- fpgen$GetFingerprint(mol[[1]])
  
}


tanimoto_similarity_df = tanimoto_distance_matrix2(fingerprints)
tanimoto_similarity_df <- as.data.frame(tanimoto_similarity_df)
#print(dim(tanimoto_similarity_df))
tanimoto_distance <- 1 - tanimoto_similarity_df

colnames(tanimoto_distance) <- merged_df$Name
rownames(tanimoto_distance) <- merged_df$Name

clusters = cluster_fingerprints(fingerprints, cutoff=0.5)




ylgnbu_col <- sequential_hcl(9, "YlGnBu")


par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)  

if (nrow(tanimoto_distance) == 1) {
  print("something")
}else{
  hmap <- Heatmap(
    as.matrix(tanimoto_distance),
    show_row_names = T,
    show_column_names = F,
    #row_names_gp = gpar(col = ifelse(rownames(distance_matrix) %in% (filter(merged_data, type == "new"))$`chemical names`, "blue", "black")),
    cluster_rows = T,
    cluster_columns = T,
    show_column_dend = T,
    show_row_dend = T,
    row_dend_reorder = T,
    column_dend_reorder = T,
    clustering_method_rows = "ward.D2", #top_annotation = colAnn,
    clustering_method_columns = "ward.D2",col = ylgnbu_col,
    width = unit(150, "mm"),border = TRUE,
    heatmap_legend_param = list(
      title = "Tanimoto Distance",
      title_position = "leftcenter-rot",
      labels_gp = gpar(fontsize = 12),
      title_gp = gpar(fontsize = 12)),
    column_gap=unit(1, "mm"))
  
  #ht = draw(hmap, heatmap_legend_side="left", annotation_legend_side="bottom")
  
  
  
}




if (nrow(tanimoto_distance) == 1) {
  heatmap_cannot_be_plotted(TRUE)
  heatmap_slot(NULL)
} else {
  hmap <- Heatmap(
    as.matrix(tanimoto_distance),
    show_row_names = TRUE,
    show_column_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_dend_reorder = TRUE,
    column_dend_reorder = TRUE,
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    col = ylgnbu_col,
    width = unit(150, "mm"),
    border = TRUE,
    heatmap_legend_param = list(
      title = "Tanimoto Distance",
      title_position = "leftcenter-rot",
      labels_gp = gpar(fontsize = 12),
      title_gp = gpar(fontsize = 12)
    ),
    column_gap = unit(1, "mm")
  )
  
  heatmap_cannot_be_plotted(FALSE)
  heatmap_slot(hmap)
}















if (nrow(distance_matrix) == 1) {
  print("display a specific message")
}else if (nrow(distance_matrix) <= 3 & nrow(distance_matrix) > 1){
  set.seed(4242)
  metaMDS_results <- metaMDS(comm = distance_matrix,
                             autotransform = FALSE,
                             engine = "monoMDS",
                             k = 2,
                             weakties = FALSE,
                             model = "global",
                             maxit = 400,
                             try = 40,
                             trymax = 100)
  
}else{
  set.seed(4242)
  metaMDS_results <- metaMDS(comm = distance_matrix,
                             autotransform = FALSE,
                             engine = "monoMDS",
                             k = 3,
                             weakties = FALSE,
                             model = "global",
                             maxit = 400,
                             try = 40,
                             trymax = 100)
  
  
}else {
  plotlyOutput("imageOutput3", width = 800, height = 800)
}


metaMDS_results.points <- data.frame(metaMDS_results$points)


metaMDS_results.points$cluster_ID <- merged_df$cluster_membership

metaMDS_results.points$Name <- merged_df$Name
#print(head(metaMDS_results.points))
metaMDS_results.points <- metaMDS_results.points %>% mutate(cluster_ID = as.character(cluster_ID))

metaMDS_results.points <- metaMDS_results.points[order(metaMDS_results.points$cluster_ID, method = "shell"),]

metaMDS_results.points$cluster_ID <- factor(metaMDS_results.points$cluster_ID,
                                            levels = as.character(1:length(clusters)), ordered = TRUE)

clustered_mol <- metaMDS_results.points %>% group_by(cluster_ID) %>% 
  summarize(cluster_count = sum(cluster_ID >= 1, na.rm = TRUE)) %>% 
  filter(cluster_count > 1)


metaMDS_results.points$cluster_ID_toPlot <- ifelse(metaMDS_results.points$cluster_ID %in% clustered_mol$cluster_ID, 
                                                   "plot", "do not plot")


metaMDS_results.points %>% 
  ggplot(., aes(x = MDS1, y = MDS2)) +
  geom_point(size = 2, col = "grey") +
  #geom_point(data = filter(metaMDS_results.points, cluster_ID_toPlot == "plot"), aes(x = MDS1, y = MDS3, col = cluster_ID), size = 4) +
  #scale_color_manual(values = c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#bf5b17','#666666')) +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  labs(title = "NMDS Plot from Distance Matrix",
       x = "MDS1",
       y = "MDS2")



















###############################################################################
### saving user input files to a folder
###############################################################################



library(shiny)
library(uuid)

ui <- fluidPage(
  fileInput("file_upload", "Upload a file"),
  textOutput("upload_status")
)

server <- function(input, output, session) {
  
  observeEvent(input$file_upload, {
    req(input$file_upload)
    
    # Generate a unique filename using UUID and timestamp
    session_id <- session$token
    original_name <- input$file_upload$name
    new_filename <- paste0("upload_", session_id, "_", original_name)
    
    # Define destination path (make sure this directory exists and is writable)
    dest_dir <- "user_uploads"
    if (!dir.exists(dest_dir)) dir.create(dest_dir)
    dest_path <- file.path(dest_dir, new_filename)
    
    # Copy the uploaded file from temp to permanent storage
    file.copy(from = input$file_upload$datapath, to = dest_path)
    
    
  })
}

shinyApp(ui, server)





















