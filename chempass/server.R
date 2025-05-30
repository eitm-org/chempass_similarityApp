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

#' rdkit conda environment
use_miniconda("my-rdkit-env2")

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
  
  #' at the start of the app, all buttons are disabled along with drop down menus
  observe({
    shinyjs::disable("process_file")
    shinyjs::disable("fingerprint_type")
    shinyjs::disable("fingerprint_button")
    shinyjs::disable("cluster")
    shinyjs::disable("cutoff")
    
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
      mixed_input <- read.csv(input$file_upload$datapath, header = FALSE, stringsAsFactors = FALSE) 
      colnames(mixed_input) <- "V1"
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
      clusters_slot(NULL)
      NMDS_slot(NULL)
      NMDS_cannot_be_plotted <- reactiveVal(FALSE)
      heatmap_cannot_be_plotted <- reactiveVal(FALSE)
    }
    
    
    
    
  })
  
  ###############################################################################
  ### processing user input data - connecting with PubChem
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
      if (ncol(mixed_input) != 1) {
        sendSweetAlert(
          session,
          title = "Incorrect format",
          text = paste("Expected 1 column, but found ", ncol(mixed_input)),
          type = "error",
          #timer = 4000,
          #showConfirmButton = FALSE
        )
        shinyjs::reset("file_upload")
        return(NULL)
      }
      
      #' trims the white spaces around text in input
      mixed_input$V1 <- trimws(mixed_input$V1)
      
      #' identifies duplicated rows
      #' and saves them to a reactive element for use later
      duplicate_rows <- mixed_input[duplicated(mixed_input), , drop = FALSE]
      duplicate_rows_slot(duplicate_rows)
      ######################
      # duplicates message in app
      ######################
      
      
      if (nrow(duplicate_rows) > 0) {
        duplicate_rows_slot(duplicate_rows)
        sendSweetAlert(
          session,
          title = "Duplicates Detected",
          text = paste("Found ", nrow(duplicate_rows), " duplicate entries in your data. Processing will continue."),
          type = "warning",
          position = "center",
          timer = 3000,
          btn_labels = NA
        )
      } else {
        sendSweetAlert(
          session,
          title = "No duplicates found",
          text = "Your data looks clean! Processing will continue.",
          type = "success",
          position = "center",
          timer = 3000,
          btn_labels = NA
        )
        #shinyjs::reset("file_upload")
      }
      
      
      mixed_input <- mixed_input %>% distinct(V1, .keep_all = TRUE)
      
      
    }, error = function(e) {
      sendSweetAlert(
        session,
        title = "Error reading file",
        text = paste("Please ensure the file is a valid plain-text CSV or TXT.\n\n"),
        type = "error",
        #timer = 6000,
        #showConfirmButton = FALSE
      )
      shinyjs::reset("file_upload")
    })
    
    shinyjs::disable("process_file")
    
    
    show_modal_spinner(text = "Processing input: Pulling data from PubChem, might take >1 min. for large lists of compounds...")
    
     
    numbers_only <- function(x) !grepl("\\D", x)
    
    mixed_input$DTXSID <- grepl("DTXSID", mixed_input$V1)
    
    mixed_input$CID <- numbers_only(mixed_input$V1)
    
    dtxsid_df <- filter(mixed_input, DTXSID == TRUE)
    dtxsid_df <- subset(dtxsid_df, select = V1)
    
    
    cid_df <- filter(mixed_input, CID == TRUE)
    cid_df <- subset(cid_df, select = V1)
    
    
    smiles_df <- filter(mixed_input, DTXSID == FALSE & CID == FALSE)
    smiles_df <- subset(smiles_df, select = V1)
    
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
    #' merged_df <- cid_processing(dfs_list)
    
    df <- bind_rows(dfs_list)
    
    
    
    #' checking for errors which are added during Pubchem API search
    error_rows <- df %>%
      filter(if_any(everything(), ~ . == 'error'))
    
    remove_modal_spinner()
    
    
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
      
      
      missing_data(c(error_rows$input, missing_names))
      
      #' is this needed if I use clean names somwhere else
      #' check the difference here
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
      
      
      
    }
    
    
    
    
    
  }) ## ending second observation here = for user provided input
  
  
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
    
    property_df <- get_properties_from_CIDs(paste(df_filt$CID, collapse = ","))
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
      title = "Please check Heatmap tab for results",
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
      
      heatmap_slot(hmap)
      heatmap_cannot_be_plotted(FALSE)
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
  ### generating clusters based on user specified cutoffs
  ###############################################################################
  
  observeEvent(input$cluster, {
    
    sendSweetAlert(
      session,
      title = "Please check Butina cluster and NMDS tab for results",
      type = "success",
      position = "center",
      timer = 3000,
      btn_labels = NA
    )
    
    req(event3_trigger())
    req(input$cutoff)
    fingerprints <- fingerprint_slot()
    merged_df <- merged_matrix()
    shinyjs::disable("cluster")
    #shinyjs::disable("cutoff")
    
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
    
    
    to_save <- merged_df[,-grep("ROMol", colnames(merged_df))]
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
                                 trymax = 100)
      
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
                                 trymax = 100)
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
  display_image <- function(file_path, title) {
    img_r <- readPNG(file_path)
    raster_img <- rasterGrob(img_r, interpolate = TRUE,
                             width = unit(0.7, "npc"),
                             height = unit(0.7, "npc"))
    
    arranged <- arrangeGrob(
      raster_img,
      top = textGrob(title, gp = gpar(fontsize = 12, col = "black"))
    )
    
    return(arranged)
  }
  
  output$cluster_pdf <- renderUI({
    
    if (is.null(clusters_slot())) {
      div(
        style = "text-align: center; padding-top: 50px;",
        h4("Butina Clusters will be generated after clustering is performed...")
      )
    } else {
      # Create a temporary file path for the PDF
      pdf_path <- tempfile(fileext = ".pdf")
      
      # Generate multipage PDF
      indices_greater_than1 <- indices_greater_than1_slot()
      
      if (length(indices_greater_than1) == 0) {
        return(
          tags$div(
            style = "padding: 20px; color: #b00;",
            tags$h4("No clusters found."),
            tags$p("Try adjusting the cutoff parameter or check your input.")
          )
        )
      }
      
      cluster_grobs <- lapply(seq_along(indices_greater_than1), function(i) {
        file_path <- paste0("cluster", i, ".png")
        title <- paste("Cluster", i)
        display_image(file_path, title)
      })
      
      pdf(pdf_path, width = 12, height = 8)
      for (grob in cluster_grobs) {
        grid.newpage()
        grid.draw(grob)
      }
      dev.off()
      
      # Move the PDF into www/ so it can be served
      www_pdf_path <- file.path("www", "clusters_only.pdf")
      file.copy(pdf_path, www_pdf_path, overwrite = TRUE)
      
      # Embed the PDF using iframe
      tags$iframe(style = "height:80vh; width:100%;", src = "clusters_only.pdf")
    }
    
  })
  
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
      
      #' creating a temporary directory
      temp_dir <- tempdir()
      
      #' saving the molecular properties
      molecular_properties <- file.path(temp_dir, "molecular-properties.csv")
      write.csv(results_data(), molecular_properties, row.names = FALSE)
      
      #' saving the second csv of failed IDs
      missing_ID <- file.path(temp_dir, "Failed_ID_from_pubchem.csv")
      write.table(missing_data(), missing_ID, row.names = FALSE, col.names = FALSE,
                quote = FALSE, sep = ",")
      
      #' Save the cluster PDF
      cluster_pdf <- file.path(temp_dir, paste0("Butina_clusters_", input$fingerprint_type, "_",
                                                input$cutoff, ".pdf"))
      
      if (input$fingerprint_type == "ECFP4") {
        save_clusters_to_pdf(clusters_slot(), mol_struct(), merged_matrix(), cluster_pdf)
      }else{
        save_clusters_to_pdf_fcfp4(clusters_slot(), mol_struct(), merged_matrix(), cluster_pdf)
      }
      
      #' saving the nmds plots
      #' 
      NMDS_png <- file.path(temp_dir, "NMDS_plot.png")
      if (NMDS_cannot_be_plotted()) {
        print("not really making an nmds plot")
      } else{
        
        png(NMDS_png, units = "in", res = 600, width = 8, height = 8)
        plot(NMDS_slot())
        dev.off()
        
      }
      
      
      
      #' Save the heatmap png
      heatmap_png <- file.path(temp_dir, "heatmap.png")
      
      if (heatmap_cannot_be_plotted()) {
        print("not really making a heatmap")
      }else{
        png(heatmap_png, units = "in", res = 600, width = 20, height = 12)
        draw(heatmap_slot(), heatmap_legend_side="left", annotation_legend_side="bottom")
        dev.off()
      }
      
      
      
      
      
      #' Providing a README of molecular properties
      readme_path <- file.path(temp_dir, "README.txt")
      writeLines(readme_text, con = readme_path)
      
      
      #' providing any duplicate rows
      duplicate_path <- file.path(temp_dir, "duplicated_entries.csv")
      write.table(duplicate_rows_slot(), duplicate_path, row.names = FALSE, col.names = FALSE,
                quote = FALSE, sep = ",")
      
      #' Providing a README of inputs from user
      if (input$use_example == TRUE) {
        inputs <- c(input$fingerprint_type, input$cutoff, input$use_example)
        names(inputs) <- c("Selected fingerprint type:", "Chosen clustering cutoff:",
                           "Used example data:")
      }else{
        
        inputs <- c(input$fingerprint_type, input$cutoff)
        names(inputs) <- c("Selected fingerprint type:", "Chosen clustering cutoff:")
      }
      
      inputs_path <- file.path(temp_dir, "User_provided_inputs.txt")
      write.table(inputs, inputs_path, sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
      
      #' creating the zip file
      zip(file, c(molecular_properties, missing_ID, cluster_pdf, heatmap_png, NMDS_png, readme_path,
                  inputs_path, duplicate_path),
          flags = "-j")
      
    },
    contentType = "text/csv/pdf/png"
  )
  

  
  
}
