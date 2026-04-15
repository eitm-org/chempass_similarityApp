# server.R

library(reticulate)
library(png)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
library(colorspace)
library(vegan)
library(base64enc)
library(readxl)
library(tibble)
library(dplyr)
library(Matrix)
library(scales)
library(stringr)
library(ggplot2)
library(plotly)
library(tidyr) ########### this needs to be installed probably in Docker??? check later

#' rdkit conda environment
#use_miniconda("my-rdkit-env2")

#' python functions that are used throughout the script
#' mainly contains functions for getting information from Pubchem
source_python("functions.py")


#' importing python libraries
pd <- import("pandas")
requests <- import("requests")
tqdm <- import("tqdm")
np <- import("numpy")
Chem <- import("rdkit.Chem")
AllChem <- import("rdkit.Chem.AllChem")
Descriptors <- import("rdkit.Chem.Descriptors")
DataStructs <- import("rdkit.DataStructs")
Butina <- import("rdkit.ML.Cluster.Butina")
PandasTools <- import("rdkit.Chem.PandasTools")
plt <- import("matplotlib.pyplot")
Draw <- import("rdkit.Chem.Draw")

#' R functions like tanimoto similarity
#' cleaning compound names
#' saving cluster images to pdf
#' readme text
source("R_functions.R")

#' start of the server 
function(input, output, session) {
  
  #' defining reactive elements here at various stages
  #' 
  reactive_df <- reactiveVal(NULL)
  event2_trigger <- reactiveVal(FALSE) 
  event3_trigger <- reactiveVal(FALSE)
  clustering_done <- reactiveVal(FALSE)
  NMD_plotting_done <- reactiveVal(FALSE)
  tanimoto_matrix <- reactiveVal(NULL)
  merged_matrix <- reactiveVal(NULL)
  results_data <- reactiveVal(NULL)
  missing_data <- reactiveVal(NULL)
  clusters_slot <- reactiveVal(NULL)
  mol_struct <- reactiveVal(NULL)
  fingerprint_slot <- reactiveVal(NULL)
  heatmap_slot <- reactiveVal(NULL)
  NMDS_slot <- reactiveVal(NULL)
  indices_greater_than1_slot <- reactiveVal(NULL)
  duplicate_rows_slot <- reactiveVal(NULL)
  toggle <- reactiveVal(NULL)
  metaMDS_results.points_slot <- reactiveVal(NULL)
  NMDS_cannot_be_plotted <- reactiveVal(FALSE)
  heatmap_cannot_be_plotted <- reactiveVal(FALSE)
  reactive_df_Ref <- reactiveVal(NULL)
  reactive_df_target <- reactiveVal(NULL)
  reactive_df_result <- reactiveVal(NULL)
  reactive_df_comb <- reactiveVal(NULL)
  event4_trigger <- reactiveVal(FALSE) 
  reactive_df_clusters <- reactiveVal(NULL)
  reactive_df_clusters_download <- reactiveVal(NULL)
  scores_df_slot <- reactiveVal(NULL)
  max_score_cutoff_slot <- reactiveVal(NULL)
  noPubchem_slot <- reactiveVal(FALSE)
  silhouette_cannot_be_plotted <- reactiveVal(FALSE)
  invalid_rows_slot <- reactiveVal(NULL)
  comb_to_save_slot <- reactiveVal(NULL)
  invalid_rows_comb_slot <- reactiveVal(NULL)
  no_pubchem_df_save_slot <- reactiveVal(NULL)
  plot1_slot <- reactiveVal(NULL)
  plot2_slot <- reactiveVal(NULL)
  plot3_slot <- reactiveVal(NULL)
  values <- reactiveValues()
  ref_duplicates_slot <- reactiveVal(NULL)
  target_duplicates_slot <- reactiveVal(NULL)
  ref_example <- reactiveVal(NULL)
  target_example <- reactiveVal(NULL)
  toggle2_list <- reactiveVal(NULL)
  #' at the start of the app, all buttons are disabled along with drop down menus
  observe({
    shinyjs::disable("process_file")
    shinyjs::disable("fingerprint_type")
    shinyjs::disable("fingerprint_button")
    shinyjs::disable("cluster")
    shinyjs::disable("cutoff")
    shinyjs::disable("process_two_files")
    shinyjs::disable("fingerprint_type_2")
    #shinyjs::disable("cluster_2")
    
    #print(input$noPubChem[1])
    
    
  })
  
  #' this observation chunk handles the example data processing

  observeEvent(input$use_example, {
    
    if(input$use_example == TRUE){
      print("user set toggle to on")
      toggle <- "on"
      mixed_input <- read.csv("CID_example_data.csv", header = FALSE, stringsAsFactors = FALSE)
      colnames(mixed_input) <- "V1"
      reactive_df(mixed_input)
      shinyjs::disable("process_file")
      shinyjs::disable("fingerprint_type")
      shinyjs::disable("fingerprint_button")
      shinyjs::disable("cluster")
      shinyjs::disable("cutoff")
      
      event2_trigger(FALSE)
      event3_trigger(FALSE)
      clustering_done(FALSE)
      NMD_plotting_done(FALSE)
      heatmap_slot(NULL)
      clusters_slot(NULL)
      NMDS_slot(NULL)
      NMDS_cannot_be_plotted <- reactiveVal(FALSE)
      heatmap_cannot_be_plotted <- reactiveVal(FALSE)
      
    }else{
      #' this resets the file upload 
      #' so if user sets toggle on and the example data sets processed
      #' and then user sets the toggle off again, it resets the whole app
      print("user set toggle to off")
      reset("file_upload")
      shinyjs::disable("process_file")
      shinyjs::disable("fingerprint_type")
      shinyjs::disable("fingerprint_button")
      shinyjs::disable("cluster")
      shinyjs::disable("cutoff")
      
      event2_trigger(FALSE)
      event3_trigger(FALSE)
      clustering_done(FALSE)
      NMD_plotting_done(FALSE)
      heatmap_slot(NULL)
      scores_df_slot(NULL)
      clusters_slot(NULL)
      NMDS_slot(NULL)
      NMDS_cannot_be_plotted <- reactiveVal(FALSE)
      heatmap_cannot_be_plotted <- reactiveVal(FALSE)
      
    }
    
    
    
  
  })
  
  
  #' this observeEvent is for checking uploaded file format
  #' files need to be either csv or txt
  #' check if someone does not save an extension to file then what happens here?
  observeEvent(input$file_upload, {
    req(input$file_upload)
    
    file_ext <- tools::file_ext(input$file_upload$name)
    print(file_ext)
    
    if (!(file_ext %in% c("csv", "txt"))) {
      sendSweetAlert(
        session,
        title = "Unsupported file type",
        text = "Please upload a .csv or .txt file only.",
        type = "error",
        #timer = 4000,
        #showConfirmButton = FALSE
      )
      shinyjs::reset("file_upload")
      return(NULL)
    }else{
      
      sendSweetAlert(
        session,
        title = "Correct file type uploaded.",
        text = tags$span(
          "Please click the", tags$b("Process This File"), "button to proceed"
        ),
        html = TRUE,
        #text = "Please click the Process This File button to proceed",
        type = "success",
        #timer = 3000,
        #showConfirmButton = FALSE
      )
      
      #' once the file format check happens the files gets read in as mixed_input
      #' still all buttons and dropdown menus are greyed out
      mixed_input <- read.csv(input$file_upload$datapath, header = T) 
      #colnames(mixed_input) <- "V1"
      reactive_df(mixed_input)
      shinyjs::enable("process_file")
      shinyjs::disable("fingerprint_type")
      shinyjs::disable("fingerprint_button")
      shinyjs::disable("cluster")
      shinyjs::disable("cutoff")
      
      event2_trigger(FALSE)
      event3_trigger(FALSE)
      clustering_done(FALSE)
      NMD_plotting_done(FALSE)
      heatmap_slot(NULL)
      scores_df_slot(NULL)
      clusters_slot(NULL)
      NMDS_slot(NULL)
      NMDS_cannot_be_plotted <- reactiveVal(FALSE)
      heatmap_cannot_be_plotted <- reactiveVal(FALSE)
    }
    
    
    
    
  })
  
  ###############################################################################
  ### processing user input data 
  ###############################################################################
  
  observeEvent(input$process_file, {
    
    req(reactive_df())
    mixed_input <- reactive_df()
    
    
    ######################
    # catching incorrect formats
    ######################
    
    tryCatch({
      
      #' Check number of columns
      #' if user provides more than 1 column than it errors out and resets file uplaod
      #' what if the other column is just NAs which sometimes happens with csv in excel sheets
      #' test this
      #' ideally the tool should be able to remove that column it its all NAs
      
      given_columns <- colnames(mixed_input)
      required_smiles <- any(tolower(given_columns) %in% c("cid", "smiles", "dtxsid", "id"))
      has_names <- any(tolower(given_columns) %in% c("name", "names"))
      valid_columns <- all(tolower(given_columns) %in% c("id","cid", "dtxsid","smile", "smiles", "name", "names"))
      
      if (!required_smiles || !valid_columns) {
        sendSweetAlert(
          session,
          title = "Incorrect columns",
          text = paste("Expected 'ID/SMILES' column and optionally 'Names'. Found:", paste(given_columns, collapse = ", ")),
          type = "error"
        )
        shinyjs::reset("file_upload")
        return(NULL)
      }
      
      # pull out blank rows
      blank_rows <- apply(mixed_input, 2, function(col) sum(is.na(col) | col == ""))
      cols_with_blanks <- names(blank_rows[blank_rows > 0])
      
      
      if (length(cols_with_blanks) > 0) {
        sendSweetAlert(
          session,
          title = "Blank rows detected",
          text = paste("Columns with blank entries:", paste(cols_with_blanks, collapse = ", ")),
          type = "error"
        )
        shinyjs::reset("file_upload")
        return(NULL)
      }
      
      # change the column names 
      colnames(mixed_input) <- case_when(
        tolower(colnames(mixed_input)) %in% c("id","cid", "smiles", "dtxsid") ~ "ID",
        tolower(colnames(mixed_input)) %in% c("name", "names") ~ "Name"
        
      )
      
      
      
      #' trims the white spaces around text in input
      mixed_input$ID <- trimws(mixed_input$ID)
      mixed_input$Name <- trimws(mixed_input$Name)
      
      
      ######################
      # duplicates message in app
      ######################
      # Check for duplicates
      duplicated_smiles <- mixed_input[duplicated(mixed_input[["ID"]]), ]
      
      
      if (NROW(duplicated_smiles) > 0) {
        
        mixed_input <- mixed_input[!duplicated(mixed_input[["ID"]]), ]
        
        sendSweetAlert(
          session,
          title = "Duplicates Detected",
          text = paste("Found duplicate IDs:", nrow(duplicated_smiles), ". Duplicate entries will be removed."),
          type = "warning",
          position = "center",
          timer = 3000,
          btn_labels = NA
        )
        
      }
      
      # Check for duplicates
      duplicated_names <- mixed_input[duplicated(mixed_input[["Name"]]), ]
      
      if (NROW(duplicated_names) > 0) {
        
        mixed_input <- mixed_input[!duplicated(mixed_input[["Name"]]), ]
        
        sendSweetAlert(
          session,
          title = "Duplicates Detected",
          text = paste("Found duplicate Names:", nrow(duplicated_names), ". Duplicate entries will be removed."),
          type = "warning",
          position = "center",
          timer = 3000,
          btn_labels = NA
        )
        
        
      } 
      
      
      duplicates <- bind_rows(duplicated_smiles, duplicated_names)
      
      
      duplicate_rows_slot(duplicates)
      #' identifies duplicated rows
      #' and saves them to a reactive element for use later
      duplicate_rows <- mixed_input[duplicated(mixed_input), , drop = FALSE]
      duplicate_rows_slot(duplicate_rows)
      
      
      
    })
    
    shinyjs::disable("process_file")
    
    
    
    
    #' this is where it takes the mixed input column after file format and duplication check
    #' and alerts the user if there are duplicates and stores them in a slot
    
    reactive_df(mixed_input)
    
    shinyjs::enable("process_file_withPubChem")
    shinyjs::enable("noPubChem")
    
  })
  
  
  observeEvent(input$process_file_withPubChem, {
    
    
    show_modal_spinner(text = "Processing input: Pulling data from PubChem, might take >1 min. for large lists of compounds...")
    
    mixed_input <- reactive_df()
    numbers_only <- function(x) !grepl("\\D", x)
    
    mixed_input$DTXSID <- grepl("DTXSID", mixed_input$ID)
    
    
    
    mixed_input$CID <- numbers_only(mixed_input$ID)
    
    dtxsid_df <- filter(mixed_input, DTXSID == TRUE)
    dtxsid_df <- subset(dtxsid_df, select = ID)
    
    print("this should print dtxsids")
    print(head(dtxsid_df))
    
    cid_df <- filter(mixed_input, CID == TRUE)
    cid_df <- subset(cid_df, select = ID)
    
    
    smiles_df <- filter(mixed_input, DTXSID == FALSE & CID == FALSE)
    smiles_df <- subset(smiles_df, select = ID)
    
    #' these individual components could be written up as functions
    #' so that they can be individually tested with testthat
    #' 
    #' This could be a dtxsid function
    #' usage eg: 
    #' df1 <- dtxsid_pubchem(dtxsid_df)
    #' 
    
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
    
    df <- bind_rows(dfs_list)
    
    
    
    #' checking for errors which are added during Pubchem API search
    error_rows <- df %>%
      filter(if_any(everything(), ~ . == 'error'))
    
    remove_modal_spinner()
    
    print("printint out error rows before heatmap plotting")
    print(error_rows)
    print(df)
    
    if (nrow(error_rows) > 0) {
      sendSweetAlert(
        session,
        title = "Errors found:",
        text = paste0("Found ", nrow(error_rows), " rows with no matches in Pubchem. Processing will continue."),
        type = "warning",
        timer = 3000,
        btn_labels = NA
      )
    }
    
    if (nrow(error_rows) == nrow(df)) {
      sendSweetAlert(
        session,
        title = "Error from Pubchem",
        text = "Input had no matched CIDs, please inspect input list",
        type = "error",
        #timer = 6000,
        #showConfirmButton = FALSE
      )
      shinyjs::reset("file_upload")
      return(NULL)
    }else{
      
      
      show_modal_spinner(text = "Processing input: Pulling data from PubChem, might take >1 min. for large lists of compounds...")
      
      #print("this process keeps running")
      
      if (dim(error_rows)[1] > 0) {
        df_filt <- filter(df, !(input %in% error_rows$input))
        
      }else{
        df_filt <- df}
      
      
      
      #' this removes the title column since it gets added during cid check
      df_filt <- subset(df_filt, select = -c(TITLE))
      
      
      #' mutates the cid into a character
      
      df_filt <- mutate(df_filt, CID = as.character(CID))
      
      df_filt <- df_filt %>% distinct(CID, .keep_all = TRUE)
      
      # get properties from Pubchem
      property_df <- tryCatch({
        get_properties_from_CIDs(df_filt$CID)
      }, error = function(e) {
        #message("Retrying once after error...")
        Sys.sleep(1)
        tryCatch(get_properties_from_CIDs(df_filt$CID), error = function(e2) {
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
      
      
      if (any(colnames(merged_df) %in% c("Name"))) {
        merged_df <- merged_df %>% select(-c(Title))
      }else{
        
        colnames(merged_df)[grep("Title", colnames(merged_df))] <- "Name"
      }
      
      
      
      
      
      #' saves the error and missing CIDs and names
      #' is it possible to show here the missing input instead
      #' might be better for user
      
      
      missing_data(c(error_rows$input, missing_names))
      
      print("missing data after pressing pubchem button")
      print(missing_data())
      
      
      #' is this needed if I use clean names somwhere else
      #' check the difference here
      merged_df$Name <- gsub(" ", "_", merged_df$Name)
      
      
      
    }
    
    
    #shinyjs::disable("noPubChem")
    
    
    
    merged_df$ROMol <- lapply(merged_df$SMILES, function(smiles) {
      
      return(Chem$MolFromSmiles(smiles))
    })
    
    merged_df <- as.data.frame(merged_df)
    
    
    
    shinyjs::disable("process_file_withPubChem")
    shinyjs::disable("noPubChem")
    
    
    merged_matrix(merged_df)
    mol_struct(merged_df$ROMol)
    event2_trigger(TRUE)
    shinyjs::enable("fingerprint_button")
    shinyjs::enable("fingerprint_type")
    
    rm(mixed_input)
    
    remove_modal_spinner()
    print(noPubchem_slot())
    
    
    
  }) ### end event for processing with PubChem
  
  
  observeEvent(input$noPubChem, {
    
    ## add in chemical formula from rdkit as a column
    # allow name column from user
    
    mixed_input <- reactive_df()
    
    print("This is after pressing the no_pubchem button blah")
    print(mixed_input)
    #colnames(mixed_input) <- "SMILES"
    
    # change the column names 
    colnames(mixed_input) <- case_when(
      tolower(colnames(mixed_input)) %in% c("id") ~ "SMILES",
      tolower(colnames(mixed_input)) %in% c("name", "names") ~ "Name"
      
    )
    
    merged_df <- mixed_input
    
    print("This is after pressing the no_pubchem button")
    print(merged_df)
    
    
    merged_df$ROMol <- lapply(merged_df$SMILES, function(smiles){
      
      return(Chem$MolFromSmiles(smiles)) 
      
    })
    
    merged_df <- as.data.frame(merged_df)
    
    #print("This is after generating structures for all input smiles")
    #print(head(merged_df))
    
    
    valid_rows <- !sapply(merged_df$ROMol, is.null)
    
    invalid_rows <- sapply(merged_df$ROMol, is.null)
    
    invalid_rows <- merged_df[invalid_rows,]
    
    
    
    if(dim(invalid_rows)[1] > 0){
      sendSweetAlert(
        session,
        title = "Error in parsing some input SMILES; these will be removed from analysis",
        type = "warning",
        position = "center",
        timer = 3000,
        btn_labels = NA
      )
      
    }
    
    invalid_rows <- invalid_rows %>% select(SMILES)
    
    ## keeping valid rows in dataframe for downstream analysis
    merged_df <- merged_df[valid_rows,]
    
    if (any(colnames(merged_df) %in% c("Name"))) {
      print("do nothing")
    }else{
      
      merged_df$Name <- paste0("SMILES_", 1:nrow(merged_df))
    }
    
    
    
    
    
    # saving the molecular properties to a csv
    mol_prop_df <- calculate_mol_mlwght_formula(merged_df)
    
    
    mol_prop_df$ROMol <- sapply(mol_prop_df$ROMol, toString)
    
    
    no_pubchem_df_save_slot(mol_prop_df)
    
    
    
    
    noPubchem <- "no_pubchem"
    
    shinyjs::disable("process_file_withPubChem")
    shinyjs::disable("noPubChem")
    
    merged_matrix(merged_df)
    mol_struct(merged_df$ROMol)
    event2_trigger(TRUE)
    invalid_rows_slot(invalid_rows)
    shinyjs::enable("fingerprint_button")
    shinyjs::enable("fingerprint_type")
    noPubchem_slot(noPubchem)
    rm(mixed_input)
    
    remove_modal_spinner()
    
  }) ## No pubchem ending second observation here = for user provided input
  
  
  ###############################################################################
  ### input example data processing - connecting with PubChem
  ###############################################################################
  
  
  observeEvent(input$use_example, {
    req(reactive_df())
    req(input$use_example)
    mixed_input <- reactive_df()
    
    mixed_input$V1 <- trimws(mixed_input$V1)
    
    
    ##duplicate_rows <- mixed_input[duplicated(mixed_input),]
    #duplicate_rows_slot(duplicate_rows)
    #mixed_input <- mixed_input %>% count(V1) %>% filter(n == 1) %>% subset(select = V1)
    
    
    show_modal_spinner(text = "Processing example CIDs: Pulling data from PubChem...")
    
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
        cid <- get_cid_from_dtxsid(dtxsid)
        
        pubchem_dtx <- append(pubchem_dtx, list(list(input = dtxsid,'DTXSID' = dtxsid, 'CID' = cid)))
        
      }
      
      df1 <- pd$DataFrame(pubchem_dtx)
    }
    
    
    
    if (dim(smiles_df)[1] > 0) {
      pubchem_smiles <- list()
      smiles <- smiles_df$V1
      
      for (smile in smiles) {
        encoded <- base64encode(charToRaw(smile))
        cid <- get_CID_from_SMILES(encoded)
        
        pubchem_smiles <- append(pubchem_smiles, list(list(input = smile,'DTXSID' = NULL, 'CID' = cid)))
        
      }
      
      df2 <- pd$DataFrame(pubchem_smiles)
    }
    
    if (dim(cid_df)[1] > 0) {
      pubchem_cid <- list()
      cids <- cid_df$V1
      
      for (cid in cids) {
        
        pubchem_cid <- append(pubchem_cid, list(list(input = cid, 'DTXSID' = NULL, 'CID' = cid)))
        
      }
      
      df3 <- pd$DataFrame(pubchem_cid)
    }
    
    
    dfs_list <- list()
    if (exists("df1")) dfs_list$df1 <- df1
    if (exists("df2")) dfs_list$df2 <- df2
    if (exists("df3")) dfs_list$df3 <- df3
    
    
    df <- bind_rows(dfs_list)
    
    error_rows <- df %>%
      filter(if_any(everything(), ~ . == 'error'))
    
    df_filt <- filter(df, !(input %in% error_rows$input))
    df_filt <- mutate(df_filt, CID = as.character(CID))
    print(str(df_filt))
    
    property_df <- get_properties_from_CIDs(df_filt$CID)
    
    
    property_df <- mutate(property_df, CID = as.character(CID))
    
    synonym_list <- list()
    for (cid in df_filt$CID) {
      synonyms <- get_synonyms_from_CID(cid)
      
      synonym_list <- append(synonym_list, list(list('CID' = cid,'SYNONYMS' = synonyms)))
      
    }
    syn_df <- pd$DataFrame(synonym_list)
    syn_df <- mutate(syn_df, CID = as.character(CID))
    
    
    df_filt <- left_join(df_filt, property_df)
    df_filt <- left_join(df_filt, syn_df)
    
    
    
    merged_df <- df_filt
    #print(merged_df)
    colnames(merged_df)[grep("Title", colnames(merged_df))] <- "Name"
    #print(merged_df)
    
    missing_data(error_rows$input)
    merged_df$Name <- gsub(" ", "_", merged_df$Name)
    merged_df$ROMol <- lapply(merged_df$SMILES, function(smiles) {
      
      return(Chem$MolFromSmiles(smiles))
    })
    
    merged_df <- as.data.frame(merged_df)
    
    #print(merged_df)
    
   
    merged_matrix(merged_df)
    mol_struct(merged_df$ROMol)
    event2_trigger(TRUE)
    shinyjs::enable("fingerprint_button")
    shinyjs::enable("fingerprint_type")
    rm(mixed_input)
    
    remove_modal_spinner()
    
  }) ## ending second observation here = for example data
  
  
  ###############################################################################
  ### generating fingerprints either ECFP/FCFP4 
  ###############################################################################
  
  observeEvent(input$fingerprint_button, {
    
    sendSweetAlert(
      session,
      title = "Please check Heatmap + Silhouette Scores tab for results",
      type = "success",
      position = "center",
      timer = 3000,
      btn_labels = NA
    )
    
    req(event2_trigger())
    req(input$fingerprint_type)
    merged_df <- merged_matrix()
    shinyjs::disable("fingerprint_button")
    #shinyjs::disable("fingerprint_type")
    
    ######################
    # fingerprint generation
    ######################
    
    
    
    if (input$fingerprint_type == "ECFP4") {
      fpgen <- AllChem$GetMorganGenerator(radius = as.integer(2), fpSize = as.integer(2048))
      fingerprints <- list()
      for (i in 1:nrow(merged_df)){ 
        mol <- merged_df$ROMol[i]
        fingerprints[[i]] <- fpgen$GetFingerprint(mol[[1]])
        
      }
      
    }else{
      invgen = AllChem$GetMorganFeatureAtomInvGen()
      ffpgen = AllChem$GetMorganGenerator(radius=as.integer(2), atomInvariantsGenerator=invgen)
      
      fingerprints <- list()
      for (i in 1:nrow(merged_df)){ 
        mol <- merged_df$ROMol[i]
        fingerprints[[i]] <- ffpgen$GetFingerprint(mol[[1]])
        
      }
    }
    
    #print(input$fingerprint_type)
    
    
    tanimoto_similarity_df = tanimoto_distance_matrix2(fingerprints)
    tanimoto_similarity_df <- as.data.frame(tanimoto_similarity_df)
    #print(dim(tanimoto_similarity_df))
    tanimoto_distance <- 1 - tanimoto_similarity_df
    
    colnames(tanimoto_distance) <- merged_df$Name
    rownames(tanimoto_distance) <- merged_df$Name
    
    tanimoto_matrix(tanimoto_distance)
    
    fingerprint_slot(fingerprints)
    
    
    ######################
    # heatmap generation
    ######################
    ylgnbu_col <- sequential_hcl(9, "YlGnBu")
    
    
    par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)  
    
    print(nrow(tanimoto_distance))
    
    
    if (nrow(tanimoto_distance) == 1) {
      heatmap_cannot_be_plotted(TRUE)
      heatmap_slot(TRUE)
    }else{
      
      hmap <- Heatmap(
        as.matrix(tanimoto_distance),
        show_row_names = ifelse(nrow(tanimoto_distance) >= 500, FALSE, TRUE),
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
      
      heatmap_slot(hmap)
      heatmap_cannot_be_plotted(FALSE)
    }
    
    
    ######################
    # silhouette scores
    ######################
    scores <- run_silhouette_analysis(fingerprints)
    scores_df <- as.data.frame(unlist(scores))
    
    print("this is after silhouette scores are generated")
    print(scores_df)
    
    if(dim(scores_df)[1] == 0){
      
      print("silhouette scores cannot be generated")
      
      sendSweetAlert(
        session,
        title = "No silhouette scores will be generated",
        type = "warning",
        position = "center",
        timer = 3000,
        btn_labels = NA
      )
      
      scores_df_slot(TRUE)
      silhouette_cannot_be_plotted(TRUE)
      
      
    }else{
      
      colnames(scores_df) <- "score"
      scores_df <- rownames_to_column(scores_df, "cutoff")
      scores_df <- mutate(scores_df, cutoff = as.numeric(cutoff))
      
      
      
      cluster_counts <- get_cluster_counts(fingerprints)
      cluster_counts <- as.data.frame(unlist(cluster_counts))
      colnames(cluster_counts) <- "cluster_counts"
      cluster_counts <- rownames_to_column(cluster_counts, "cutoff")
      cluster_counts <- mutate(cluster_counts, cutoff = as.numeric(cutoff))
      
      
      scores_df <- left_join(scores_df, cluster_counts)
      
      
      
      max_score_cutoff <- scores_df %>% slice_max(order_by = score, n = 1) %>% 
        pull(cutoff)
      
      max_score <- round(max(scores_df$score), 2)
      
      plot1 <- ggplot(scores_df, aes(x = cutoff, y = score)) +
        geom_point(size = 4) +
        geom_line() +
        geom_vline(xintercept = max_score_cutoff, color = "red", linetype = "dashed",
                   size = 1.2) +
        scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.05)) +
        theme_bw() +
        theme(text = element_text(size = 15),
              panel.grid = element_line(colour = "grey")) +
        labs(x = "Clustering Cutoff (Distance)", y = "Silhouette Score") +
        ggtitle(paste("Max Silhouette Score = ", max_score, "; at distance cutoff = ", max_score_cutoff))
      
      plot2 <- ggplot(scores_df, aes(x = cutoff, y = cluster_counts)) +
        geom_point(size = 4) +
        geom_line() +
        scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.05)) +
        theme_bw() +
        theme(text = element_text(size = 15),
              panel.grid = element_line(colour = "grey")) +
        labs(x = "Clustering Cutoff (Distance)", y = "Number of Clusters")
      
      plot3 <- ggplot(scores_df, aes(x = cluster_counts, y = score)) +
        geom_point(aes(col = cutoff), size = 4) +
        theme_bw() +
        theme(text = element_text(size = 15),
              panel.grid = element_line(colour = "grey")) +
        scale_colour_viridis_c(option = "viridis") +
        labs(y = "Silhouette Score", x = "Number of Clusters")
      
      
      plot1_slot(plot1)
      plot2_slot(plot2)
      plot3_slot(plot3)
      
      scores_df_slot(TRUE)
      silhouette_cannot_be_plotted(FALSE)
      
    }
    
    
    
    
    
    event3_trigger(TRUE)
    shinyjs::enable("fingerprint_button")
    shinyjs::enable("cutoff")
    shinyjs::enable("cluster")
    NMDS_slot(NULL)
    clusters_slot(NULL)
    
    
  }) # ending the third observation event for fingerprint type
  
  ###############################################################################
  ### Rendering the heatmap based on the fingerprints and distance matrix generated above
  ###############################################################################
  
  
  output$imageOutput1 <- renderPlot({
    
    
    
    req(heatmap_slot())
    draw(heatmap_slot(), heatmap_legend_side="left", annotation_legend_side="bottom")
    
    
    
  }, res = 100)
  
  
  output$heatmap_container <- renderUI({
    if (is.null(heatmap_slot())) {
      div(
        style = "text-align: center; padding-top: 50px;",
        h4("Heatmap will be generated after Fingerprint generation...")
      )
    } else if(isTRUE(heatmap_cannot_be_plotted())){
      div(
        style = "text-align: center; padding-top: 50px;",
        h4("Heatmap cannot be plotted with a single data point...")
      )
      
    }else {
      plotOutput("imageOutput1", width = 1500, height = 1500)
    }
  })
  
  
  ###############################################################################
  ### Rendering the silhouette scores
  ###############################################################################
  
  
  output$imageOutput_SScore <- renderPlot({
    
    req(plot1_slot())
    req(plot2_slot())
    req(plot3_slot())
    
    gridExtra::grid.arrange(plot1_slot(),plot2_slot(),plot3_slot(),ncol = 1)
    
  })
  
  
  output$silhouette_container <- renderUI({
    if (is.null(scores_df_slot())) {
      div(
        style = "text-align: center; padding-top: 50px;",
        h4("Silhouette scores will be generated after Fingerprint generation...")
      )
    } else if(isTRUE(silhouette_cannot_be_plotted())){
      div(
        style = "text-align: center; padding-top: 50px;",
        h4("No clusters for silhouette score plotting...")
      )
      
    } 
    else {
      plotOutput("imageOutput_SScore", width = 1500, height = 1500)
    }
  })
  
  
  
  
  ###############################################################################
  ### generating clusters based on user specified cutoffs
  ###############################################################################
  
  observeEvent(input$cluster, {
    
    
    
    req(event3_trigger())
    req(input$cutoff)
    fingerprints <- fingerprint_slot()
    merged_df <- merged_matrix()
    shinyjs::disable("cluster")
    #shinyjs::disable("cutoff")
    
    show_modal_spinner(text = "Generating cluster pdf...\nPlease Download results and plots")
    clusters = cluster_fingerprints(fingerprints, cutoff=input$cutoff)
    
    clusters_slot(clusters)
    
    
    merged_df$cluster_membership <- 1000
    
    for (i in seq_along(clusters)) {
      cluster_indices <- clusters[[i]]
      cluster_mols <- merged_df$Name[unlist(cluster_indices) + 1]
      merged_df$cluster_membership <- ifelse(merged_df$Name %in% cluster_mols, i, merged_df$cluster_membership)
    }
    
    ######################
    # saving cluster png
    ######################
    
    indices_greater_than1 <- which(sapply(clusters, length) > 1)
    clustered_mols <- list()
    images <- list()
    for (i in seq_along(indices_greater_than1)) {
      
      cluster_indices <- clusters[[i]]
      
      clustered_mols[[i]] <- merged_df$ROMol[unlist(cluster_indices) + 1]
      opts = Draw$MolDrawOptions()
      opts$legendFraction = 0.2
      opts$legendFontSize = 20L
      img <- Draw$MolsToGridImage(clustered_mols[[i]], subImgSize = list(400L, 400L),
                                  drawOptions=opts,
                                  legends=merged_df$Name[unlist(cluster_indices) +1])
      img$save(paste0("cluster", i, ".png"))
      
      
      images[[i]] <- img
    }
    
    indices_greater_than1_slot(indices_greater_than1)
    merged_matrix(merged_df)
    
    to_save <- merged_df
    #to_save <- merged_df[,-grep("ROMol", colnames(merged_df))]
    to_save$ROMol <- sapply(merged_df$ROMol, toString)
    results_data(to_save)
    
    clustering_done(TRUE)
    shinyjs::enable("cluster")
    
    distance_matrix <- tanimoto_matrix()
    
    ######################
    # NMDS analysis
    ######################
    
    if (nrow(distance_matrix) == 1) {
      NMDS_cannot_be_plotted(TRUE)
      NMDS_slot(TRUE)
      NMD_plotting_done(TRUE)
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
                                 trymax = 50)
      
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
        filter(cluster_count > 1) %>% ungroup()
      
      
      metaMDS_results.points$cluster_ID_toPlot <- ifelse(metaMDS_results.points$cluster_ID %in% clustered_mol$cluster_ID, 
                                                         "plot", "do not plot")
      
      
      
      ######################
      # NMDS plotting- ggplot
      ######################
      
      plot <- metaMDS_results.points %>% 
        ggplot(., aes(x = MDS1, y = MDS2)) +
        geom_point(size = 2, col = "grey") +
        geom_point(data = filter(metaMDS_results.points, cluster_ID_toPlot == "plot"), aes(x = MDS1, y = MDS2, col = cluster_ID), size = 4) +
        #scale_color_manual(values = c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#bf5b17','#666666')) +
        theme_classic() +
        theme(text = element_text(size = 15)) +
        labs(title = "NMDS Plot from Distance Matrix",
             x = "MDS1",
             y = "MDS2")
      
      NMDS_slot(plot)
      
      NMD_plotting_done(TRUE)
      metaMDS_results.points_slot(metaMDS_results.points)
      
      
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
                                 trymax = 50)
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
        filter(cluster_count > 1) %>% ungroup()
      
      
      metaMDS_results.points$cluster_ID_toPlot <- ifelse(metaMDS_results.points$cluster_ID %in% clustered_mol$cluster_ID, 
                                                         "plot", "do not plot")
      
      
      
      ######################
      # NMDS plotting- ggplot
      ######################
      
      plot <- metaMDS_results.points %>% 
        ggplot(., aes(x = MDS1, y = MDS2)) +
        geom_point(size = 2, col = "grey") +
        geom_point(data = filter(metaMDS_results.points, cluster_ID_toPlot == "plot"), aes(x = MDS1, y = MDS2, col = cluster_ID), size = 4) +
        #scale_color_manual(values = c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#bf5b17','#666666')) +
        theme_classic() +
        theme(text = element_text(size = 15)) +
        labs(title = "NMDS Plot from Distance Matrix",
             x = "MDS1",
             y = "MDS2")
      
      NMDS_slot(plot)
      
      NMD_plotting_done(TRUE)
      metaMDS_results.points_slot(metaMDS_results.points)
      
    }
    
    ###############
    ## cluster pdf rendering
    ###############
    temp_dir1 <- tempfile()
    dir.create(temp_dir1)
    pdf_path <- file.path(temp_dir1, "clusters_only.pdf")
    
    
    if (input$fingerprint_type == "ECFP4") {
      save_clusters_to_pdf(clusters_slot(), mol_struct(), merged_matrix(), pdf_path)
      
    }else{
      save_clusters_to_pdf_fcfp4(clusters_slot(), mol_struct(), merged_matrix(), pdf_path)
      
    }
    
    
    
    values$cluster_pdf_path <- pdf_path
    values$singleList_pdf <- TRUE
    
    addResourcePath(paste0("pdfs_", session$token), temp_dir1)
    
    remove_modal_spinner()
    
    
    sendSweetAlert(
      session,
      title = "Please check Butina cluster and NMDS tab for results",
      type = "success",
      position = "center",
      timer = 3000,
      btn_labels = NA
    )
    
  }) # ending the fourth observation of cluster cutoff
  
  ###############################################################################
  ### Rendering the NMDS plots - interactive based on distance matrix and clustering cutoffs
  ###############################################################################
  
  output$imageOutput3 <- renderPlotly({
    
    req(event3_trigger())
    req(tanimoto_matrix())
    req(clustering_done())
    req(NMDS_slot())
    
    metaMDS_results.points <- metaMDS_results.points_slot()
    
    
    ######################
    # NMDS - interactive plots
    ######################
    
    g <- metaMDS_results.points %>% 
      ggplot(., aes(x = MDS1, y = MDS2, text = Name)) +
      geom_point(size = 2, col = "grey") +
      geom_point(data = filter(metaMDS_results.points, cluster_ID_toPlot == "plot"), aes(x = MDS1, y = MDS2, col = cluster_ID), size = 4) +
      #scale_color_manual(values = c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#bf5b17','#666666')) +
      theme_classic() +
      theme(text = element_text(size = 15)) +
      labs(title = "NMDS Plot from Distance Matrix",
           x = "MDS1",
           y = "MDS2")
    
    ggplotly(g, tooltip = "text")
    
    
  })
  
  output$NMDS_container <- renderUI({
    if (is.null(NMDS_slot())) {
      div(
        style = "text-align: center; padding-top: 50px;",
        h4("NMDS plot will be generated after clustering is performed...")
      )
    }else if (NMDS_cannot_be_plotted()) {
      div(
        style = "text-align: center; padding-top: 50px;",
        h4("NMDS plot cannot be generated, data points N = 1...")
      )
      
    }else {
      plotlyOutput("imageOutput3", width = 800, height = 800)
    } 
    
    
  })
  
  ###############################################################################
  ### Rendering the Butina clusters - needs Butina clustering to be done
  ###############################################################################
  
  
  ######################
  # function - image display
  ######################
  output$cluster_pdf <- renderUI({
    
    if (is.null(clusters_slot())) {
      div(
        style = "text-align: center; padding-top: 50px;",
        h4("Butina Clusters will be generated after clustering is performed...")
      )
    } else {
      
      req(values$singleList_pdf)
      
      # Embed the PDF using iframe
      tags$iframe(style = "height:80vh; width:100%;", src = paste0("pdfs_", session$token, "/clusters_only.pdf"))
    }
    
  })
  
  ###############################################################################
  ### the two list comparison pdf as a separate tab
  ###############################################################################
  
  
  
  
  ### allowing user to upload example file
  
  observeEvent(input$example_btn, {
    req(input$example_btn)
    
    # Path to your example CSV
    example1_path <- "Reference_isomers_enentiomers.csv"  
    example2_path <- "Target_isomers_enentiomers.csv"  
    
    
    
    if (file.exists(example1_path)) {
      reactive_df_Ref(read.csv(example1_path))
      print(head(reactive_df_Ref()))
      
    } else {
      showNotification("Reference file not found", type = "error")
    }
    if (file.exists(example2_path)) {
      reactive_df_target(read.csv(example2_path))
      print(head(reactive_df_target()))
    } else {
      showNotification("Target file not found", type = "error")
    }
    
    
    
    shinyjs::reset("reference_file_upload")
    shinyjs::reset("target_file_upload")
    shinyjs::disable("example_btn")
    shinyjs::enable("process_two_files")
    
  }) ## end of observation event for checking loading example reference and target files
  
  
  
  ### generating the observe event ; waiting for input files from user
  
  #observeevent here later
  
  observeEvent(input$reference_file_upload, {
    
    #req(input$reference_file_upload)
    
    file_ext_Ref <- tools::file_ext(input$reference_file_upload$name)
    if (!(file_ext_Ref %in% c("csv", "txt"))) {
      sendSweetAlert(
        session,
        title = "Unsupported reference file type",
        text = "Please upload a .csv or .txt file only.",
        type = "error"
        #timer = 4000,
        #showConfirmButton = FALSE
      )
      shinyjs::reset("reference_file_upload")
      return(NULL)
    }else{
      
      #' once the file format check happens the files gets read in as mixed_input
      #' still all buttons and dropdown menus are greyed out
      ref_list <- read.csv(input$reference_file_upload$datapath, stringsAsFactors = FALSE) 
      
      # pull out blank rows
      blank_rows <- apply(ref_list, 2, function(col) sum(is.na(col) | col == ""))
      cols_with_blanks <- names(blank_rows[blank_rows > 0])
      
      
      if (length(cols_with_blanks) > 0) {
        sendSweetAlert(
          session,
          title = "Blank rows detected",
          text = paste("Columns with blank entries:", paste(cols_with_blanks, collapse = ", ")),
          type = "error"
        )
        shinyjs::reset("reference_file_upload")
        return(NULL)
      }
      
      ## handling column name issues
      given_columns <- colnames(ref_list)
      required_smiles <- any(tolower(given_columns) %in% c("smile", "smiles"))
      has_names <- any(tolower(given_columns) %in% c("name", "names"))
      valid_columns <- all(tolower(given_columns) %in% c("smile", "smiles", "name", "names"))
      
      
      if (!required_smiles || !valid_columns) {
        sendSweetAlert(
          session,
          title = "Incorrect columns",
          text = paste("Expected 'SMILES' column and optionally 'Names'. Found:", paste(given_columns, collapse = ", ")),
          type = "error"
        )
        shinyjs::reset("reference_file_upload")
        return(NULL)
      }
      
      if (has_names) {
        sendSweetAlert(
          session,
          title = "Reference SMILES and compound Names provided by user",
          type = "success"
        )
      }
      
      
      # change the column names 
      colnames(ref_list) <- case_when(
        tolower(colnames(ref_list)) %in% c("smile", "smiles") ~ "SMILES",
        tolower(colnames(ref_list)) %in% c("name", "names") ~ "Name"
      
      )
      
      
      # Check for duplicates
      duplicated_smiles <- ref_list[duplicated(ref_list[["SMILES"]]), ]
      
      
      if (NROW(duplicated_smiles) > 0) {
        
        ref_list <- ref_list[!duplicated(ref_list[["SMILES"]]), ]
        
        sendSweetAlert(
          session,
          title = "Duplicates Detected",
          text = paste("Found duplicate SMILES:", nrow(duplicated_smiles), ". Duplicate entries will be removed. Please upload target list."),
          type = "warning",
          position = "center",
          timer = 3000,
          btn_labels = NA
        )
        
      }
      
      # Check for duplicates
      duplicated_names <- ref_list[duplicated(ref_list[["Name"]]), ]
      
      if (NROW(duplicated_names) > 0) {
        
        ref_list <- ref_list[!duplicated(ref_list[["Name"]]), ]
        
        sendSweetAlert(
          session,
          title = "Duplicates Detected",
          text = paste("Found duplicate Names:", nrow(duplicated_names), ". Duplicate entries will be removed. Please upload target list."),
          type = "warning",
          position = "center",
          timer = 3000,
          btn_labels = NA
        )
        
        
      } 
      
      
      duplicates <- bind_rows(duplicated_smiles, duplicated_names)
      
      ref_list$group <- "reference"
      reactive_df_Ref(ref_list) 
      ref_duplicates_slot(duplicates)
      shinyjs::disable("example_btn")
      
    }
    
    
    
  })
  
  
  
  
  observeEvent(input$target_file_upload, {
    
    req(input$target_file_upload)
    
    file_ext_target <- tools::file_ext(input$target_file_upload$name)
    if (!(file_ext_target %in% c("csv", "txt"))) {
      sendSweetAlert(
        session,
        title = "Unsupported target file type",
        text = "Please upload a .csv or .txt file only.",
        type = "error",
        #timer = 4000,
        #showConfirmButton = FALSE
      )
      shinyjs::reset("target_file_upload")
      return(NULL)
    }else{
      
      #' once the file format check happens the files gets read in as mixed_input
      #' still all buttons and dropdown menus are greyed out
      target_list <- read.csv(input$target_file_upload$datapath, stringsAsFactors = FALSE) 
      
      # pull out blank rows
      blank_rows <- apply(target_list, 2, function(col) sum(is.na(col) | col == ""))
      cols_with_blanks <- names(blank_rows[blank_rows > 0])
      
      
      if (length(cols_with_blanks) > 0) {
        sendSweetAlert(
          session,
          title = "Blank rows detected",
          text = paste("Columns with blank entries:", paste(cols_with_blanks, collapse = ", ")),
          type = "error"
        )
        shinyjs::reset("target_file_upload")
        return(NULL)
      }
      
      
      
      ## handling column name issues
      given_columns <- colnames(target_list)
      required_smiles <- any(tolower(given_columns) %in% c("smile", "smiles"))
      has_names <- any(tolower(given_columns) %in% c("name", "names"))
      valid_columns <- all(tolower(given_columns) %in% c("smile", "smiles", "name", "names"))
      
      if (!required_smiles || !valid_columns) {
        sendSweetAlert(
          session,
          title = "Incorrect columns",
          text = paste("Expected 'SMILES' column and optionally 'Names'. Found:", paste(given_columns, collapse = ", ")),
          type = "error"
        )
        shinyjs::reset("target_file_upload")
        return(NULL)
      }
      
      if (has_names) {
        sendSweetAlert(
          session,
          title = "Target SMILES and compound Names provided by user",
          type = "success"
        )
      }
      
      # change the column names 
      colnames(target_list) <- case_when(
        tolower(colnames(target_list)) %in% c("smile", "smiles") ~ "SMILES",
        tolower(colnames(target_list)) %in% c("name", "names") ~ "Name"
        
      )
      
      
      # check for SMILES duplicates
      duplicated_smiles <- target_list[duplicated(target_list[["SMILES"]]), ]
      
      # Check for duplicates
      if (NROW(duplicated_smiles) > 0) {
        
        target_list <- target_list[!duplicated(target_list[["SMILES"]]), ]
        
        sendSweetAlert(
          session,
          title = "Duplicates Detected",
          text = paste("Found duplicate SMILES:", nrow(duplicated_smiles), ". Duplicate entries will be removed. Please upload target list."),
          type = "warning",
          position = "center",
          #timer = 3000,
          #btn_labels = NA
        )
        
      }
      
      # check for Names duplicates
      duplicated_names <- target_list[duplicated(target_list[["Name"]]), ]
      
      if (NROW(duplicated_names) > 0) {
        
        target_list <- target_list[!duplicated(target_list[["Name"]]), ]
        
        sendSweetAlert(
          session,
          title = "Duplicates Detected",
          text = paste("Found duplicate Names:", nrow(duplicated_names), ". Duplicate entries will be removed. Please upload target list."),
          type = "warning",
          position = "center",
          #timer = 3000,
          #btn_labels = NA
        )
        
        
      } 
      
      duplicates <- bind_rows(duplicated_smiles, duplicated_names)
      
      target_list$group <- "target"
      reactive_df_target(target_list)
      target_duplicates_slot(duplicates)
      shinyjs::enable("process_two_files")
      shinyjs::disable("example_btn")
    }
    
    
    
  }) ## end of observation event for checking file type; target list
  
  
  
  observeEvent(input$process_two_files, {
    req(reactive_df_Ref())
    req(reactive_df_target())
    
    shinyjs::disable("process_two_files")
  
    
    
    ref_df <- reactive_df_Ref()
    target_df <- reactive_df_target()
    comb <- rbind(ref_df, target_df)
    
    
    duplicate_smiles <- (comb %>% count(SMILES) %>% filter(n > 1))$SMILES
    
    comb <- comb %>% distinct(SMILES, .keep_all = TRUE)
    comb <- comb %>% distinct(Name, .keep_all = TRUE)
  
    # from SMILES; compare lists
    
    comb$ROMol <- lapply(comb$SMILES, function(smiles) {
      
      return(Chem$MolFromSmiles(smiles))
    })
    
    comb <- as.data.frame(comb)
    
    invalid_rows_comb <- sapply(comb$ROMol, is.null)
    
    invalid_rows_comb <- comb[invalid_rows_comb,]
    
    
    
    if(NROW(invalid_rows_comb) > 0){
      sendSweetAlert(
        session,
        title = "Error in parsing some input SMILES; these will be removed from analysis",
        type = "warning",
        position = "center",
        timer = 3000,
        btn_labels = NA
      )
      
    }
    
    invalid_rows_comb <- invalid_rows_comb %>% select(SMILES)
    
    valid_rows_comb <- !sapply(comb$ROMol, is.null)
    comb <- comb[valid_rows_comb,]
    
    print(head(comb))
    
    length_ref <- length((filter(comb, group == "reference"))$SMILES)
    length_target <- length((filter(comb, group == "target"))$SMILES)
    
    
    
    if (any(colnames(comb) %in% c("Name"))) {
      print("do nothing")
    }else{
      
      comb$Name <- "random"
      
      comb$Name[grep("reference", comb$group)] <- paste0("Ref_SMILES_", 1:length_ref)
      comb$Name[grep("target", comb$group)] <- paste0("Target_SMILES_", 1:length_target)
      
      
    }
    
   
    
    
    comb_mol_prop <- calculate_mol_mlwght_formula(comb)
    
    
    
    reactive_df_comb(comb_mol_prop)
    invalid_rows_comb_slot(invalid_rows_comb)
    event4_trigger(TRUE)
    shinyjs::enable("fingerprint_button_2")
    shinyjs::enable("fingerprint_type_2")
     
  }) ## end of event for combining and processing input files
  
  observeEvent(input$fingerprint_button_2,{
    
    req(reactive_df_comb())
    req(input$fingerprint_button_2)
    req(input$fingerprint_type_2)
    
    comb <- reactive_df_comb()
    
    ref_list_name <- (comb %>% filter(group == "reference"))$Name
    target_list_name <- (comb %>% filter(group == "target"))$Name
    
    shinyjs::disable("fingerprint_button_2")
    #############################################
    ## fingerprint conversion; give the user a choice
    #############################################
    
    if (input$fingerprint_type_2 == "ECFP4") {
      fpgen <- AllChem$GetMorganGenerator(radius = as.integer(2), fpSize = as.integer(2048))
      fingerprints_ref_target <- list()
      for (i in 1:nrow(comb)){ 
        mol <- comb$ROMol[i]
        fingerprints_ref_target[[i]] <- fpgen$GetFingerprint(mol[[1]])
        
      }
      
    }else{
      invgen = AllChem$GetMorganFeatureAtomInvGen()
      ffpgen = AllChem$GetMorganGenerator(radius=as.integer(2), atomInvariantsGenerator=invgen)
      
      fingerprints_ref_target <- list()
      for (i in 1:nrow(comb)){ 
        mol <- comb$ROMol[i]
        fingerprints_ref_target[[i]] <- ffpgen$GetFingerprint(mol[[1]])
        
      }
    }
    
    print("fingerprint generation done")
    tanimoto_similarity_df = tanimoto_distance_matrix2(fingerprints_ref_target)
    tanimoto_similarity_df <- as.data.frame(tanimoto_similarity_df)
    print(dim(tanimoto_similarity_df))
    
    print(head(tanimoto_similarity_df))
    
    colnames(tanimoto_similarity_df) <- comb$Name
    rownames(tanimoto_similarity_df) <- comb$Name
    
    # Convert similarity matrix to a tidy data frame
    similarity_long <- as.data.frame(tanimoto_similarity_df) %>%
      mutate(target = rownames(tanimoto_similarity_df)) %>%
      pivot_longer(
        cols = -target,
        names_to = "reference",
        values_to = "similarity"
      )
    
    print(head(similarity_long))
    
    #print(similarity_long)
    result <- similarity_long %>%
      filter(target %in% target_list_name, reference %in% ref_list_name) %>%
      group_by(target) %>%
      slice_max(similarity, n = 3) %>% 
      mutate(rank = row_number()) %>% 
      ungroup()
    
    print(head(result))
    
  
    sendSweetAlert(
      session,
      title = "Please check Two List Compare tab for results",
      type = "success",
      position = "center",
      timer = 3000,
      btn_labels = NA
    )
    
    reactive_df_result(result)
    
    top_hit <- result %>% group_by(target) %>% 
      filter(rank == 1) %>% ungroup()
    
    other_hits <- result %>% mutate(similarity = round(similarity, digits = 2)) %>% 
      group_by(target) %>% 
      filter(rank > 1) %>%
      summarise(other_hits = paste(reference, collapse = "; "),
                other_similarities = paste(similarity, collapse = "; ")) %>% 
      ungroup()
    
    target_df <- as.data.frame(top_hit$target)
    colnames(target_df) <- "Name"
    
    
    target_df <- left_join(target_df, comb)
    
    
    ref_df <- as.data.frame(top_hit$reference)
    colnames(ref_df) <- "Name"
    
    ref_df <- left_join(ref_df, comb)
    colnames(ref_df) <- paste0("Reference_",colnames(ref_df))
    
    
    
    similarity_df <- as.data.frame(top_hit$similarity)
    colnames(similarity_df) <- "similarity"
    
    comb_to_save <- cbind(target_df, ref_df, similarity_df) %>% 
      left_join(other_hits, by = c("Name" = "target"))
    
    
    comb_to_save <- comb_to_save[,-c(grep("group", colnames(comb_to_save)))]
    
    comb_to_save$ROMol <- sapply(comb_to_save$ROMol, toString)
    comb_to_save$Reference_ROMol <- sapply(comb_to_save$Reference_ROMol, toString)
    
  
    
    print("trying to combine the properites")
    print(head(comb_to_save))
    #print(str(comb_to_save))
    
    comb_to_save_slot(comb_to_save)
    
    #############
    # pdf rendering
    #############
    opts <- Draw$MolDrawOptions()
    opts$legendFraction <- 0.2  # Increase space for legend
    opts$legendFontSize <- 30L  # Adjust font size as needed
    
    
    temp_dir2 <- tempfile()
    dir.create(temp_dir2)
    pdf_filepath <- file.path(temp_dir2, "target_vs_reference_similarity.pdf")
    
    #pdf_path2 <- tempfile(fileext = ".pdf")
    
    pdf(pdf_filepath, width = 10, height = 8)
    
    unique_targets <- unique(result$target)
    
    for (target_name in unique_targets) {
      # Get all matches for this target (top 3)
      target_matches <- result %>%
        filter(target == target_name) %>%
        arrange(desc(similarity))
      
      # Get all reference names for this target
      reference_names <- target_matches$reference
      
      # Get similarity scores for subtitle
      scores <- round(target_matches$similarity, 3)
      
      # New plotting page
      plot.new()
      title(main = paste0("Target Molecule: ", target_name))
      
      # Filter mols: target + all reference compounds
      mols_to_plot <- filter(comb, Name %in% c(target_name, reference_names))
      
      # Reorder to ensure target is first, then references in rank order
      mols_to_plot <- mols_to_plot %>%
        mutate(order = case_when(
          Name == target_name ~ 0,
          Name == reference_names[1] ~ 1,
          Name == reference_names[2] ~ 2,
          Name == reference_names[3] ~ 3,
          TRUE ~ 4
        )) %>%
        arrange(order)
      
      legends <- sapply(1:nrow(mols_to_plot), function(i) {
        if (mols_to_plot$Name[i] == target_name) {
          paste0(mols_to_plot$Name[i], " (Target)")
        } else {
          rank <- which(reference_names == mols_to_plot$Name[i])
          score <- scores[rank]
          paste0(mols_to_plot$Name[i], " (", score, ")")
        }
      })
      
      # Create image (adjust molsPerRow if needed - 2 or 4 depending on layout preference)
      img <- Draw$MolsToGridImage(
        mols_to_plot$ROMol,
        molsPerRow = 2L,  # or 4L to show all in one row
        subImgSize = list(600L, 600L),  # adjust size as needed
        drawOptions = opts,
        legends = legends
      )
      
      # Save to temp file and draw
      temp_file <- tempfile(fileext = ".png")
      img$save(temp_file)
      plot_img <- png::readPNG(temp_file)
      graphics::rasterImage(plot_img, 0, 0.1, 1, 1)
    }
    
    dev.off()
    addResourcePath(paste0("pdfs_", session$token), temp_dir2)
    
    values$pdf_filepath <- pdf_filepath
    values$pdf_ready <- TRUE
    
    print("done making pdf of top 3 hits")
    
    shinyjs::enable("fingerprint_button_2")
  }) ### end of event after pressing fingerprint generation button
  
  
  ###############################################################################
  ### Rendering the pdf for target reference list comparison
  ###############################################################################
  
  output$two_list_pdf <- renderUI({
    
    if (is.null(reactive_df_result())) {
      div(
        style = "text-align: center; padding-top: 50px;",
        h4("PDF will be created once reference and target lists are uploaded (second tab)...")
      )
    }else{
      req(values$pdf_ready)  # Wait until PDF is generated
      
      # Embed the PDF using iframe
      tags$iframe(style = "height:80vh; width:100%;", src = paste0("pdfs_", session$token, "/target_vs_reference_similarity.pdf"))
      
    }
     
  })
  
  output$download_ui2 <- renderUI({
    req(reactive_df_result())
    
    downloadButton("downloadData2", "Download Zip File", style = "color: white; background-color: #006d2c; border-color: black;")
  }) ## adding in a download button once clustering is done
  
  output$downloadData2 <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".zip", sep = "")
    },
    
    content = function(file) {
      
      
      # creating a file list
      file_list2 <- c()
      
      #' creating a temporary directory
      temp_dir2 <- tempdir()
      
      #' saving the initial combined list from two uploaded files
      comb_list <- file.path(temp_dir2, "molecular_properties_clusterMembership.csv")
      write.csv(comb_to_save_slot(), comb_list, row.names = FALSE)
      file_list2 <- c(file_list2, comb_list)
      
      # saving the invalid SMILES
      
      if(NROW(invalid_rows_comb_slot()) == 0) {
        
        data_to_write_comb <- tibble(NA)
        
      }else{
        
        data_to_write_comb <- invalid_rows_comb_slot()
      }
      
      smiles_parsing_error <- file.path(temp_dir2, "error_SMILES_parsing.csv")
      write.table(data_to_write_comb, smiles_parsing_error, row.names = F, col.names = F, quote = F, sep = ",")
      file_list2 <- c(file_list2, smiles_parsing_error)
      
      ## duplicates
      # target_duplicates_slot
      
      if(NROW(target_duplicates_slot()) == 0) {
        
        data_to_write_target <- tibble(NA)
        
      }else{
        
        data_to_write_target <- target_duplicates_slot()
      }
      
      target_duplicates <- file.path(temp_dir2, "targetList_duplicates.csv")
      write.table(data_to_write_target, target_duplicates, row.names = F, col.names = F, quote = F, sep = ",")
      file_list2 <- c(file_list2, target_duplicates)
      
      # ref_duplicates_slot
      
      if(NROW(ref_duplicates_slot()) == 0) {
        
        data_to_write_ref <- tibble(NA)
        
      }else{
        
        data_to_write_ref <- ref_duplicates_slot()
      }
      
      ref_duplicates <- file.path(temp_dir2, "ReferenceList_duplicates.csv")
      write.table(data_to_write_ref, ref_duplicates, row.names = F, col.names = F, quote = F, sep = ",")
      file_list2 <- c(file_list2, ref_duplicates)
      
      
      ## adding a README
      readme_path2 <- file.path(temp_dir2, "README.txt")
      writeLines(readme_text_twoList, con = readme_path2)
      file_list2 <- c(file_list2, readme_path2)
      
      ## saving pdf of similarity
      
      file_list2 <- c(file_list2, values$pdf_filepath)
      
      
      
      

      
      
      #' creating the zip file
      zip(file, file_list2,
          flags = "-j")
      
      
      
    },
    contentType = "text/csv/pdf/png"
    
    
    
    
    
  )
  
  
  ###############################################################################
  ### Downloading data 
  ###############################################################################
  
  output$download_ui <- renderUI({
    req(clustering_done())
    req(NMD_plotting_done())
    
    
    downloadButton("downloadData", "Download Zip File", style = "color: white; background-color: #006d2c; border-color: black;")
  }) ## adding in a download button once clustering is done
  
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      
      # creating a file list of outputs
      files_to_zip <- c()
      
      #' creating a temporary directory
      temp_dir <- tempdir()
      
      ####################################
      #' saving the molecular properties
      #' or saving cluster membership for non pubchem data
      #################################### 
      
      # i want to check if noPubchem_slot is set to a string and if it is,
      # then combine two dataframes into one and give as output
      
      if (isTRUE(noPubchem_slot() == "no_pubchem")) {
        
        
        combined_df <- left_join(results_data(), no_pubchem_df_save_slot())
        
        molecular_properties <- file.path(temp_dir, "molecular_properties_clusterMembership.csv")
        write.csv(combined_df, molecular_properties, row.names = FALSE)
        
        files_to_zip <- c(files_to_zip, molecular_properties)
        
      }else{
        ## calling it molecular properties here to maintain download code; changing name of file
        molecular_properties <- file.path(temp_dir, "molecular_properties_clusterMembership.csv")
        write.csv(results_data(), molecular_properties, row.names = FALSE)
        
        files_to_zip <- c(files_to_zip, molecular_properties)
      }
      
      
      ####################################
      ## missing data files to csv
      #################################### 
      
      
      if (isTRUE(noPubchem_slot() == "no_pubchem") && NROW(invalid_rows_slot()) > 0) {
        data_to_write <- invalid_rows_slot()
        
        
      } else if(isFALSE(noPubchem_slot() == "no_pubchem") && NROW(missing_data()) > 0){
        data_to_write <- missing_data()
        
        
      }else{
        
        data_to_write <- tibble(NA)
      }
      
      missing_ID <- file.path(temp_dir, "Failed.csv")
      write.table(data_to_write, missing_ID, row.names = FALSE, col.names = FALSE,
                  quote = FALSE, sep = ",")
      files_to_zip <- c(files_to_zip, missing_ID)
      
      
      
      ####################################
      #' Save the cluster PDF
      ####################################
      
      
      files_to_zip <- c(files_to_zip, values$cluster_pdf_path)
     
      
      ####################################
      #' saving the nmds plots
      ####################################
      
      #' 
      
      if (isTRUE(NMDS_cannot_be_plotted())) {
        print("not really making an nmds plot")
      } else{
        
        NMDS_png <- file.path(temp_dir, "NMDS_plot.png")
        png(NMDS_png, units = "in", res = 600, width = 8, height = 8)
        plot(NMDS_slot())
        dev.off()
        
        files_to_zip <- c(files_to_zip, NMDS_png)
      }
      
      
      ####################################
      #' Save the heatmap png
      ####################################
      
      
      
      if (isTRUE(heatmap_cannot_be_plotted())) {
        print("not really making a heatmap")
      }else{
        
        heatmap_png <- file.path(temp_dir, "heatmap.png")
        png(heatmap_png, units = "in", res = 600, width = 20, height = 12)
        draw(heatmap_slot(), heatmap_legend_side="left", annotation_legend_side="bottom")
        dev.off()
        
        files_to_zip <- c(files_to_zip, heatmap_png)
      }
      
      
      ####################################
      # saving the silhouette scores pngs 
      ####################################
      
      
      if (isTRUE(silhouette_cannot_be_plotted())) {
        
        print("random")
        
      }else{
        
        silhouette_png <- file.path(temp_dir, "silhouette_scores.png")
        png(silhouette_png, units = "in", res = 600, width = 20, height = 12)
        plot1_obj <- plot1_slot()
        plot2_obj <- plot2_slot() 
        plot3_obj <- plot3_slot()
        
        gridExtra::grid.arrange(plot1_obj, plot2_obj, plot3_obj, ncol = 1)
        dev.off()
        
        files_to_zip <- c(files_to_zip, silhouette_png)
        
      }
      
      ####################################
      #' Providing a README 
      ####################################
      
      readme_path <- file.path(temp_dir, "README.txt")
      writeLines(readme_text, con = readme_path)
      files_to_zip <- c(files_to_zip, readme_path)
      
      
      ####################################
      #' providing any duplicate rows
      ####################################
      
      
      print(NROW(duplicate_rows_slot()))
      
      if (NROW(duplicate_rows_slot()) > 0) {
        
        dupl_to_write <- duplicate_rows_slot()
        
      } else{
        
        dupl_to_write <- tibble(NA)
        
      }
      
      duplicate_path <- file.path(temp_dir, "duplicated_entries.csv")
      write.table(dupl_to_write, duplicate_path, row.names = FALSE, col.names = FALSE,
                  quote = FALSE, sep = ",")
      
      files_to_zip <- c(files_to_zip, duplicate_path)
      
      
      ####################################
      #' Providing user inputs in download
      ####################################
      
      
      if (input$use_example == TRUE) {
        inputs <- c(input$fingerprint_type, input$cutoff, input$use_example)
        names(inputs) <- c("Selected fingerprint type:", "Chosen clustering cutoff:",
                           "Used example data:")
      }else if (isTRUE(noPubchem_slot() == "no_pubchem")){
        
        inputs <- c("",input$fingerprint_type, input$cutoff)
        names(inputs) <- c("User bypassed PubChem; molecular properties calculated using RDKit",
                           "Selected fingerprint type:", "Chosen clustering cutoff:")
      } else{
        
        inputs <- c("",input$fingerprint_type, input$cutoff)
        names(inputs) <- c("User chose PubChem; molecular properties from PubChem",
                           "Selected fingerprint type:", "Chosen clustering cutoff:")
        
      }
      
      inputs_path <- file.path(temp_dir, "User_provided_inputs.txt")
      write.table(inputs, inputs_path, sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
      files_to_zip <- c(files_to_zip, inputs_path)
      
      ####################################
      #' creating the zip file
      ####################################
      
      #' 
      
      if (length(files_to_zip) > 0) {
        zip(file, files_to_zip, flags = "-j")
      }
      
    },
    contentType = "text/csv/pdf/png"
  )
  
  
  
  

  
  
}
