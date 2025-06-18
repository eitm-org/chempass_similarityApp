# CHEMPASS; A Tool for Chemical Structure Similarity Scoring and Clustering.


R shiny app to get physical properties of compounds from PubChem and group compounds by chemical similarity.


This repository contains R scripts, python functions and Dockerfile for deploying the app; 

Repository Structure Overview
This section gives a quick explanation of the main folders and files in this repository:

	- /.github/workflows
	This contains a .yaml file that defines how the application is deployed based on a trigger (in our case pushing changes to main branch).

	- /chempass
	This folder contains code to run an R shiny app (scripts like server.R, ui.R). It also defines the functions in python and R used by server.R script.
	/www folder contains images that are displayed on the website.
	/tests contains R scripts to test individual functions used in the app.

	- Dockerfile
	Dockerfile that defines the commands to build a Docker image with all necessary dependencies for running the R Shiny application.
	It uses the packages outlined in my-rdkit-env2.yml to build a conda environment for running RDKit and its dependancies.


CHEMPASS can be accessed at: https://chempass.emilabs.org/

How to use the tool:
-------------------------------
Upload your list of compounds in .csv/.txt format either as CID/DTXSIDs/SMILES, select your choice of fingerprint generation, set the clustering threshold, and explore cluster-level statistics and heatmaps.

zip file for download includes a .pdf of clusters with similarity scores, molecular properties of compounds, high resolution images of heatmap and NMDS plots.

Description of tool outputs:

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
  - The totality of representations are collapsed into a 2048 bit space creating the final hash for that molecule (or “fingerprint”).

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





