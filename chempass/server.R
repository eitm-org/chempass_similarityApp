# server.R

library(reticulate)
library(png)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
#library(tidyverse)
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


use_miniconda("my-rdkit-env2")

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

source("R_functions.R")

function(input, output, session) {
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
  
  observe({
    shinyjs::disable("process_file")
    shinyjs::disable("fingerprint_type")
    shinyjs::disable("fingerprint_button")
    shinyjs::disable("cluster")
    shinyjs::disable("cutoff")
    
  })
  
  reset_counter <- reactiveVal(0)

  observeEvent(input$use_example, {
    
    if(input$use_example == TRUE){
      print("user set toggle to on")
      toggle <- "on"
      mixed_input <- read.csv("CID_example_data.csv", header = FALSE)
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
      
    }else{
      print("user set toggle to off")
      reset("form_upload")
      shinyjs::disable("process_file")
      shinyjs::disable("fingerprint_type")
      shinyjs::disable("fingerprint_button")
      shinyjs::disable("cluster")
      shinyjs::disable("cutoff")
      
      event2_trigger(FALSE)
      event3_trigger(FALSE)
      clustering_done(FALSE)
      NMD_plotting_done(FALSE)
      reset('file_upload')
    }
    
    
    
  
  })
  
  # Render fileInput, refresh it when reset_counter changes
  output$file_input_ui <- renderUI({
    reset_counter()  # dependency here
    fileInput("file_upload", "Upload your file")
  })
  
  
  
  observeEvent(input$file_upload, {
    req(input$file_upload)
    mixed_input <- read.csv(input$file_upload$datapath, header = FALSE) 
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
    
  })
  
  observeEvent(input$process_file, {
    req(reactive_df())
    mixed_input <- reactive_df()
    
    
    shinyjs::disable("process_file")
    
    
    
    mixed_input <- as.data.frame(trimws(mixed_input$V1))
    colnames(mixed_input) <- "V1"
    
    duplicate_rows <- mixed_input[duplicated(mixed_input),]
    duplicate_rows_slot(duplicate_rows)
    mixed_input <- mixed_input %>% count(V1) %>% filter(n == 1) %>% subset(select = V1)
    
    
    
    
    withProgress(message = "Processing:", value = 0, {
      
      
      
      
      numbers_only <- function(x) !grepl("\\D", x)
      
      mixed_input$DTXSID <- grepl("DTXSID", mixed_input$V1)
      
      mixed_input$CID <- numbers_only(mixed_input$V1)
      
      dtxsid_df <- filter(mixed_input, DTXSID == TRUE)
      dtxsid_df <- subset(dtxsid_df, select = V1)
      
      
      cid_df <- filter(mixed_input, CID == TRUE)
      cid_df <- subset(cid_df, select = V1)
      
      
      smiles_df <- filter(mixed_input, DTXSID == FALSE & CID == FALSE)
      smiles_df <- subset(smiles_df, select = V1)
      
      incProgress(0.5, detail = "Pulling data from PubChem...")
      
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
    })
    
    
    
  }) ## ending second observation here = for user provided input
  
  
  observeEvent(input$use_example, {
    req(reactive_df())
    req(input$use_example)
    mixed_input <- reactive_df()
    
    mixed_input <- as.data.frame(trimws(mixed_input$V1))
    colnames(mixed_input) <- "V1"
    
    duplicate_rows <- mixed_input[duplicated(mixed_input),]
    duplicate_rows_slot(duplicate_rows)
    mixed_input <- mixed_input %>% count(V1) %>% filter(n == 1) %>% subset(select = V1)
    
    
    
    
    withProgress(message = "Processing:", value = 0, {
      
      
      
      
      numbers_only <- function(x) !grepl("\\D", x)
      
      mixed_input$DTXSID <- grepl("DTXSID", mixed_input$V1)
      
      mixed_input$CID <- numbers_only(mixed_input$V1)
      
      dtxsid_df <- filter(mixed_input, DTXSID == TRUE)
      dtxsid_df <- subset(dtxsid_df, select = V1)
      
      
      cid_df <- filter(mixed_input, CID == TRUE)
      cid_df <- subset(cid_df, select = V1)
      
      
      smiles_df <- filter(mixed_input, DTXSID == FALSE & CID == FALSE)
      smiles_df <- subset(smiles_df, select = V1)
      
      incProgress(0.5, detail = "Pulling data from PubChem...")
      
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
    })
    
    
    
  }) ## ending second observation here = for example data
  
  
  
  
  observeEvent(input$fingerprint_button, {
    req(event2_trigger())
    req(input$fingerprint_type)
    merged_df <- merged_matrix()
    shinyjs::disable("fingerprint_button")
    #shinyjs::disable("fingerprint_type")
    
    
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
    
    
    #print(merged_df$Name)
    
    colnames(tanimoto_distance) <- merged_df$Name
    rownames(tanimoto_distance) <- merged_df$Name
    
    #print("this is where the rowname error was")
    
    
    
    tanimoto_matrix(tanimoto_distance)
    
    
    
    fingerprint_slot(fingerprints)
    event3_trigger(TRUE)
    shinyjs::enable("fingerprint_button")
    shinyjs::enable("cutoff")
    shinyjs::enable("cluster")
  }) # ending the third observation event for fingerprint type
  
  
  
  observeEvent(input$cluster, {
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
    #reset('file_upload')
  }) # ending the fourth observation of cluster cutoff
  
  
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
      write.csv(missing_data(), missing_ID, row.names = FALSE, col.names = FALSE)
      
      #' Save the cluster PDF
      cluster_pdf <- file.path(temp_dir, paste0("Butina_clusters_", input$fingerprint_type, "_",
                                                input$cutoff, ".pdf"))
      
      if (input$fingerprint_type == "ECFP4") {
        save_clusters_to_pdf(clusters_slot(), mol_struct(), merged_matrix(), cluster_pdf)
      }else{
        save_clusters_to_pdf_fcfp4(clusters_slot(), mol_struct(), merged_matrix(), cluster_pdf)
      }
      
      #' saving the nmds plots
      NMDS_png <- file.path(temp_dir, "NMDS_plot.png")
      png(NMDS_png, units = "in", res = 300, width = 8, height = 8)
      plot(NMDS_slot())
      dev.off()
      
      #' Save the heatmap png
      heatmap_png <- file.path(temp_dir, "heatmap.png")
      png(heatmap_png, units = "in", res = 300, width = 20, height = 12)
      draw(heatmap_slot(), heatmap_legend_side="left", annotation_legend_side="bottom")
      dev.off()
      
      #' Providing a README of molecular properties
      readme_path <- file.path(temp_dir, "README.txt")
      writeLines(readme_text, con = readme_path)
      
      
      #' providing any duplicate rows
      duplicate_path <- file.path(temp_dir, "duplicated_entries.csv")
      write.csv(duplicate_rows_slot(), duplicate_path, row.names = FALSE, col.names = FALSE)
      
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
  
  
  
  
  display_image <- function(file_path, title) {
    img_r <- readPNG(file_path)
    raster_img <- rasterGrob(img_r, interpolate = TRUE,
                             width = unit(0.7, "npc"),
                             height = unit(0.7, "npc"))
    
    arranged <- arrangeGrob(
      raster_img,
      top = textGrob(title, gp = gpar(fontsize = 10, col = "black"))
    )
    
    return(arranged)
  }
  
  output$cluster_pdf <- renderUI({
    req(clustering_done())

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
    
    pdf(pdf_path, width = 8, height = 8)
    for (grob in cluster_grobs) {
      grid.newpage()
      grid.draw(grob)
    }
    dev.off()
    
    # Move the PDF into www/ so it can be served
    www_pdf_path <- file.path("www", "clusters_only.pdf")
    file.copy(pdf_path, www_pdf_path, overwrite = TRUE)
    
    # Embed the PDF using iframe
    tags$iframe(style = "height:600px; width:100%;", src = "clusters_only.pdf")
  })
  
  
  output$imageOutput3 <- renderPlot({
    req(event2_trigger())
    req(tanimoto_matrix())
    req(clustering_done())
    distance_matrix <- tanimoto_matrix()
    merged_data <- merged_matrix()
    clusters <- clusters_slot()
    
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
    
    
    metaMDS_results.points$cluster_ID <- merged_data$cluster_membership
    
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
    
    NMDS_plot <- metaMDS_results.points %>% 
      ggplot(., aes(x = MDS1, y = MDS3)) +
      geom_point(size = 2, col = "grey") +
      geom_point(data = filter(metaMDS_results.points, cluster_ID_toPlot == "plot"), aes(x = MDS1, y = MDS3, col = cluster_ID), size = 4) +
      #scale_color_manual(values = c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#bf5b17','#666666')) +
      theme_classic() +
      theme(text = element_text(size = 15)) +
      labs(title = "NMDS Plot from Distance Matrix",
           x = "MDS1",
           y = "MDS2")
    
    plot(NMDS_plot)
    
    NMDS_slot(NMDS_plot)
    
    NMD_plotting_done(TRUE)
    
  }, res = 100)
  
  
  #outputOptions(output, "imageOutput3", suspendWhenHidden = FALSE)
  
  
  
  output$imageOutput1 <- renderPlot({
    req(tanimoto_matrix()) 
    req(event3_trigger())
    distance_matrix <- tanimoto_matrix()
    merged_data <- merged_matrix()
    
    library(ComplexHeatmap)
    
    library(colorspace)
    ylgnbu_col <- sequential_hcl(9, "YlGnBu")
    
    
    par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)  
    
    hmap <- Heatmap(
      as.matrix(distance_matrix),
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
    
    ht = draw(hmap, heatmap_legend_side="left", annotation_legend_side="bottom")
    
    heatmap_slot(hmap)
    
  }, res = 100)
  
  

  
  
}
