# R_functions.R

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

#' function to calculate the Tanimoto similarity score
tanimoto_distance_matrix2 <- function(fingerprints) {
  n <- length(fingerprints)
  distance_matrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:i) {
      distance_matrix[i, j] <- DataStructs$TanimotoSimilarity(fingerprints[[i]], fingerprints[[j]])
      distance_matrix[j, i] <- distance_matrix[i, j]
    }
  }
  return(distance_matrix)
}

#' function to clean compound names, to allow R to escape special characters
sanitize_names <- function(x) {
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  x <- gsub("_+", "_", x)           
  x <- gsub("^_|_$", "", x)         
  return(x)
}


save_clusters_to_pdf <- function(clusters, list_mol, merged_df, output_filename = "clusters.pdf") {
  
  grDevices::pdf(output_filename, width = 12, height = 8)
  
  #' Iterate through the clusters
  #' This step depends on clusters generated from fingerprints and user provided cutoffs
  for (i in seq_along(clusters)) {
    cluster_indices <- clusters[[i]]
    
    #' Skip empty clusters
    if (length(cluster_indices) == 0) {
      next
    }
    
    #' Get molecules for this cluster (adjust index by +1 for R) since the cluster indices are generated with python
    #' where numbering starts with 0
    cluster_mols <- list_mol[unlist(cluster_indices) + 1]
    
    #' Configure drawing options
    opts <- Draw$MolDrawOptions()
    opts$legendFraction <- 0.2  # Increase space for legend
    opts$legendFontSize <- 40L  # Adjust font size as needed
    #opts$aspectRatio <- 1.0
    #opts$coordScale <- 1.0
    
    #' Special handling for single molecule in cluster
    #' Had to do this because it was cutting off text for singletons
    if (length(cluster_mols) == 1) {
      
      #' Get chemical names for these molecules
      chemical_names <- merged_df$Name[unlist(cluster_indices) + 1]
      #print(paste0("processing:", chemical_names))
      #' Start plotting image
      plot.new()
      title(main = paste("Cluster", i))
      
      
      img <- Draw$MolToImage(cluster_mols[[1]], 
                             subImgSize = list(600L, 600L),
                             legend=chemical_names[1],
                             drawOptions=opts)
      
      
      temp_file <- tempfile(fileext = ".png")
      img$save(temp_file)
      plot_img <- png::readPNG(temp_file)
      graphics::rasterImage(plot_img, 0.1, 0.1, 0.5, 0.5)
      file.remove(temp_file)
    } else {
      #' Create the grid image for multiple molecules
      #' So this should capture all molecules that are clustered based on user provided cutoff
      merged_df$Name <- sanitize_names(merged_df$Name)
      
      cluster_ind_fingerprints <- list()
      fpgen <- AllChem$GetMorganGenerator(radius = as.integer(2), fpSize = as.integer(2048))
      for (j in 1:length(cluster_indices)){ 
        cluster_ind_fingerprints[[j]] <- fpgen$GetFingerprint(cluster_mols[[j]])
        
      }
      
      tanimoto_similarity_df = tanimoto_distance_matrix2(cluster_ind_fingerprints)
      tanimoto_similarity_df <- as.data.frame(tanimoto_similarity_df)
      tanimoto_similarity_df <- round(tanimoto_similarity_df, 2)
      
      
      colnames(tanimoto_similarity_df) <- merged_df$Name[unlist(cluster_indices) + 1]
      tanimoto_similarity_df$Name <- merged_df$Name[unlist(cluster_indices) + 1]
      
      cluster_center <- merged_df$Name[cluster_indices[[1]] + 1]
      
      
      
      temp <- tanimoto_similarity_df %>%
        select(contains("Name"), all_of(cluster_center))
      
      temp <-temp %>%
        arrange(desc(.data[[cluster_center]]))
      
      temp <- left_join(temp, merged_df)
      
      
      #' Get chemical names for these molecules and add the similarity scores to cluster center on a new line
      
      
      final_name <- paste(temp$Name, temp[[cluster_center]], sep = " \n")
      
      #print(temp)
      #print(cluster_center)
      #print(final_name)
      
      img <- Draw$MolsToGridImage(temp$ROMol, 
                                  subImgSize = list(600L, 600L),
                                  drawOptions = opts,
                                  legends = final_name)
      
      
      temp_file <- tempfile(fileext = ".png")
      img$save(temp_file)
      plot_img <- png::readPNG(temp_file)
      
      
      plot.new()
      title(main = paste("Cluster", i))
      graphics::rasterImage(plot_img, 0, 0.1, 1, 1)
      file.remove(temp_file)
    }
  }
  
  
  grDevices::dev.off()
  
  return(output_filename)
}

#' same function as above but uses fcfp fingerprint generation
#' Could be combined with the previous function 
save_clusters_to_pdf_fcfp4 <- function(clusters, list_mol, merged_df, output_filename = "clusters.pdf") {
  
  grDevices::pdf(output_filename, width = 12, height = 8)
  
  
  for (i in seq_along(clusters)) {
    cluster_indices <- clusters[[i]]
    
    #' Skip empty clusters
    if (length(cluster_indices) == 0) {
      next
    }
    
    
    cluster_mols <- list_mol[unlist(cluster_indices) + 1]
    
    
    opts <- Draw$MolDrawOptions()
    opts$legendFraction <- 0.2  
    opts$legendFontSize <- 40L  
    
    
    if (length(cluster_mols) == 1) {
      
      
      chemical_names <- merged_df$Name[unlist(cluster_indices) + 1]
      
      plot.new()
      title(main = paste("Cluster", i))
      
      
      img <- Draw$MolToImage(cluster_mols[[1]], 
                             subImgSize = list(600L, 600L),
                             legend=chemical_names[1],
                             drawOptions=opts)
      
      
      temp_file <- tempfile(fileext = ".png")
      img$save(temp_file)
      plot_img <- png::readPNG(temp_file)
      graphics::rasterImage(plot_img, 0.1, 0.1, 0.5, 0.5)
      file.remove(temp_file)
    } else {
      
      merged_df$Name <- sanitize_names(merged_df$Name)
      
      cluster_ind_fingerprints <- list()
      invgen = AllChem$GetMorganFeatureAtomInvGen()
      ffpgen = AllChem$GetMorganGenerator(radius=as.integer(2), atomInvariantsGenerator=invgen)
      
      for (j in 1:length(cluster_indices)){ 
        cluster_ind_fingerprints[[j]] <- ffpgen$GetFingerprint(cluster_mols[[j]])
        
      }
      
      tanimoto_similarity_df = tanimoto_distance_matrix2(cluster_ind_fingerprints)
      tanimoto_similarity_df <- as.data.frame(tanimoto_similarity_df)
      tanimoto_similarity_df <- round(tanimoto_similarity_df, 2)
      
      
      colnames(tanimoto_similarity_df) <- merged_df$Name[unlist(cluster_indices) + 1]
      tanimoto_similarity_df$Name <- merged_df$Name[unlist(cluster_indices) + 1]
      
      cluster_center <- merged_df$Name[cluster_indices[[1]] + 1]
      
      temp <- tanimoto_similarity_df %>%
        select(contains("Name"), all_of(cluster_center))
      
      temp <-temp %>%
        arrange(desc(.data[[cluster_center]]))
      temp <- left_join(temp, merged_df)
      
      
      
      final_name <- paste(temp$Name, temp[[cluster_center]], sep = " \n")
      
      img <- Draw$MolsToGridImage(temp$ROMol, 
                                  subImgSize = list(600L, 600L),
                                  drawOptions = opts,
                                  legends = final_name)
      
      
      temp_file <- tempfile(fileext = ".png")
      img$save(temp_file)
      plot_img <- png::readPNG(temp_file)
      
      
      plot.new()
      title(main = paste("Cluster", i))
      graphics::rasterImage(plot_img, 0, 0.1, 1, 1)
      file.remove(temp_file)
    }
  }
  
  
  grDevices::dev.off()
  
  return(output_filename)
}


#' function for processing input DTXSID through Pubchem

dtxsid_pubchem <- function(dtxsid_df){
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
    
  }
  return(df1)
}

#' function for processing input SMILES through Pubchem
smiles_pubchem <- function(smiles_df){
  
  if (dim(smiles_df)[1] > 0) {
    pubchem_smiles <- list()
    smiles <- smiles_df$V1
    
    for (smile in smiles) {
      smile <- gsub("\\s+", "", smile)
      encoded <- base64encode(charToRaw(smile))
      #cid <- get_CID_from_SMILES(encoded)
      
      cid <- tryCatch({
        get_CID_from_SMILES(encoded)
      }, error = function(e) {
        #message("Retrying once after error...")
        Sys.sleep(1)
        tryCatch(get_CID_from_SMILES(encoded), error = function(e2) {
          #message("Second attempt failed.")
          "error"
        })
      })
      
      pubchem_smiles <- append(pubchem_smiles, list(list(input = smile,'DTXSID' = NULL, 
                                                         'CID' = cid, 'TITLE' = NULL)))
      
    }
    
    df2 <- pd$DataFrame(pubchem_smiles)
    
  }
  
  return(df2)
  
}

#' function for processing input CID 
#' this function takes the user provided CIDs and checks if they have existing titles in Pubchem
#' since the names of compounds are used downstream, this is important to be present

input_cid_reformating <- function(cid_df){
  
  if (dim(cid_df)[1] > 0) {
    pubchem_cid <- list()
    cids <- cid_df$V1
    
    for (cid in cids) {
      
      title <- tryCatch({
        check_cid(cid)
      }, error = function(e) {
        #message("Retrying once after error...")
        Sys.sleep(1)
        tryCatch(check_cid(cid), error = function(e2) {
          #message("Second attempt failed.")
          "error"
        })
      })
      
      pubchem_cid <- append(pubchem_cid, list(list(input = cid, 'DTXSID' = NULL, 'CID' = cid,
                                                   'TITLE' = title)))
      
    }
    
    df3 <- pd$DataFrame(pubchem_cid)
    
  }
  
  return(df3)
}

#' function to get synonyms from CIDs

cid_to_synonym <- function(df_filt){
  
  synonym_list <- list()
  for (cid in df_filt$CID) {
    
    synonyms <- tryCatch({
      get_synonyms_from_CID(cid)
    }, error = function(e) {
      #message("Retrying once after error...")
      Sys.sleep(1)
      tryCatch(get_synonyms_from_CID(cid), error = function(e2) {
        #message("Second attempt failed.")
        "error"
      })
    })
    
    synonym_list <- append(synonym_list, list(list('CID' = cid,'SYNONYMS' = synonyms)))
    
  }
  
  syn_df <- pd$DataFrame(synonym_list)
  syn_df <- mutate(syn_df, CID = as.character(CID))
  
  return(syn_df)
  
}

#' function for downstream analysis of cid
#' includes getting synonyms, properties from pubchem
#' handling missing titles
#' 
cid_processing <- function(dfs_list){
  
  df <- bind_rows(dfs_list)
  
  #' checking for errors which are added during Pubchem API search
  error_rows <- df %>%
    filter(if_any(everything(), ~ . == 'error'))
  
  #' this removes the errors from the data frame
  #' if there are no error rows, it just saves the same dataframe as df_filt
  
  if (dim(error_rows)[1] > 1) {
    df_filt <- filter(df, !(input %in% error_rows$input))
    
  }else{
    df_filt <- df
  }
  
  
  
  #' now to catch the missing names from titles
  #' so this checks whether a CID has a title missing
  #' and collects those missing into a new variable called missing_names
  missing_names <- (df_filt[which(is.na(df_filt$Title)),])$CID
  
  #' if there are missing_names, it removes the entire row from the datafrmae
  #' this ensures that downstream, the CIDs that get processed are clean
  
  if(length(missing_names) > 0) {
    df_filt <- df_filt[-c(which(is.na(df_filt$Title))),]
  }
  
  #' this removes the title column since it gets added during cid check
  df_filt <- subset(df_filt, select = -c(TITLE))
  
  #' mutates the cid into a character
  df_filt <- mutate(df_filt, CID = as.character(CID))
  
  
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
  
  return(df_filt)
}











#' adding a README that gets downloaded with other files 
#' explanations for the plots and files generated as outputs for this tool
readme_text <- "
Molecular Properties Explained
-------------------------------

1. XLogP (Partition Coefficient)
   - Computationally generated octanol-water partition coefficient or distribution coefficient. 
   - XLogP is used as a measure of hydrophilicity or hydrophobicity of a molecule.

2. Molecular Weight (Mol Weight)
   - The molecular weight is the sum of all atomic weights of the constituent atoms in a compound, measured in g/mol. 
   - In the absence of explicit isotope labelling, averaged natural abundance is assumed. 
   - If an atom bears an explicit isotope label, 100% isotopic purity is assumed at this location.

3. HBondDonorCount
   - Number of hydrogen-bond donors in the structure.

4. HBondAcceptorCount
   - Number of hydrogen-bond acceptors in the structure.

5. Topological Polar Surface Area (TPSA)
   - Topological polar surface area, computed by the algorithm described in the paper by Ertl et al.

These properties are critical for assessing the 'drug-likeness' and physicochemical behavior of molecules.

Fingerprints Explained
-------------------------------

How are fingerprints generated from molecular structures?
  - Each atom in a compound is iteratively accessed at a diameter of 4 neighbors to produce a unique binary representation representing that structure. 
  - The totality of representations are collapsed into a 2048 bit space creating the final hash for that molecule (or “fingerprint”)

Butina clustering Explained
-------------------------------
  - A cluster centroid is the molecule within a given cluster which has the largest number of neighbors.
  - For each molecule, all molecules with a Tanimoto distance below a given threshold are counted.
  - Molecules are sorted by their number of neighbors in descending order, so that potential cluster centroids (i.e. the compounds with the largest number of neighbors) are placed at the top of the file.
  - Starting with the first molecule (centroid) in the sorted list.
      - All molecules with a Tanimoto index above or equal to the cut-off value used for clustering then become members of that cluster (in case of similarity).
      - Each molecule that has been identified as a member of the given cluster is flagged and removed from further comparisons. Thus, flagged molecules cannot become either another cluster centroid or a member of another cluster. This process is like putting an exclusion sphere around the newly formed cluster.
      - Once the first compound in the list has found all its neighbors, the first available (i.e. not flagged) compound at the top of the list becomes the new cluster centroid.
  - The same process is repeated for all other unflagged molecules down the list.
      - Molecules that have not been flagged by the end of the clustering process become singletons.
      - Note that some molecules assigned as singletons can have neighbors at the given Tanimoto similarity index, but those neighbors have been excluded by a stronger cluster centroid.

Butina_clusters pdf Explained
-------------------------------
  - Each cluster is represented on a single page of the pdf
      - For each cluster generated with Butina clustering, intra-cluster Tanimoto similarity scores are calculated.
      - On each page, the first molecule displayed is the cluster center and the Tanimoto similarity between that cluster center and other cluster members are displayed below the compound name.
      - Typically, the minimum intra-cluster similarity score displayed here is = (1-clustering cutoff set by user).
      - Singletons are plotted per page.
      
NMDS Explained
-------------------------------
  - Non-metric MultiDimensional Scaling (NMDS) is a distance-based ordination technique.  Because it focuses on the distance matrix, it is very flexible – any distance measure can be used.
  - NMDS does not plot the distances per se among sample units.  
  - Instead, it attempts to plot sample units in such a way that the distances among sample units in the ordination space are in the same rank order as the distances among sample units as measured by the original distance matrix.
  - It does not seek to maximize the variability associated with individual axes of the ordination. As a result, the axes of an NMDS ordination are entirely arbitrary.
  
Heatmap Explained
-------------------------------
  - Heatmap is generated with Tanimoto distance measure = (1 - Tanimoto similarity score)
  - Heatmap is generated using the ComplexHeatmap() package, with clustering method for rows and columns set to ward.D2. 

"