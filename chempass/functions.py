from rdkit.Chem import AllChem
import pandas as pd
from tqdm import tqdm
from rdkit.Chem import AllChem, Descriptors
from rdkit import DataStructs
import seaborn as sn
from rdkit.ML.Cluster import Butina
import numpy as np
from rdkit.Chem import PandasTools
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import matplotlib.pyplot as plt
from rdkit.Chem import Draw
import requests
import re
from io import StringIO
import base64
import time


from sklearn.metrics import silhouette_score, silhouette_samples
#import seaborn as sns
#from collections import defaultdict


### calculating molecular properties
# calculate the logP using RDkit and add to dataframe
# function definition
def calculate_mol_mlwght_formula(data, mols_column = 'ROMol'):

    '''
    This function takes in a dataframe of Molecules; currently the column name has to be ROMol
    Then using rdkit's molecule, it calculates the exact mass and formula and adds it to another dataframe with exact_mass and Formula as column names
    This needs to happen for both the pubchem avoiding step as well as compare two lists
    Also provides Lipinski's rule of 5: XLogP, HBondDonorCount, TPSA, HBondAcceptorCount
    '''

    data_copy = data.copy()
    

    
    # empty list for calculated exact molecular weights
    exactmol_weights = []

    # empty list for calculated molecular formula
    Mol_formula = []

    # empty list for calculated XlogP
    XLog_P = []

    # empty list for calculated HBondDonorCount
    HBond_DonorCount = []

    # empty list for calculated HBondAcceptorCount
    HBond_AcceptorCount = []

    # empty list for calculated TPSA
    TPSAList = []



    for mol in data_copy[mols_column]:
        
        # Calculate molecular weight values and append to list
        ExactMolWt = rdMolDescriptors.CalcExactMolWt(mol)
        exactmol_weights.append(ExactMolWt)


        # calculate molecular formula from structure
        MolFormula = rdMolDescriptors.CalcMolFormula(mol)
        Mol_formula.append(MolFormula)

        # calculate XLogP from structure
        XlogP = Descriptors.MolLogP(mol)
        XLog_P.append(XlogP)

        # calculate HBondDonorCount from structure
        hbd_count = Descriptors.NumHDonors(mol)
        HBond_DonorCount.append(hbd_count)

        # calculate HBondAcceptorCount from structure
        hba_count = Descriptors.NumHAcceptors(mol)
        HBond_AcceptorCount.append(hba_count)

        # calculate TPSA from structure
        TPSA = Descriptors.TPSA(mol)
        TPSAList.append(TPSA)

        

    data_copy['exact_mass'] = exactmol_weights
    data_copy['Formula'] = Mol_formula
    data_copy['XLogP'] = XLog_P
    data_copy['HBondDonorCount'] = HBond_DonorCount
    data_copy['HBondAcceptorCount'] = HBond_AcceptorCount
    data_copy['TPSA'] = TPSAList


    return data_copy



def tanimoto_distance_matrix(fp_list):
    """Calculate distance matrix for fingerprint list"""
    dissimilarity_matrix = []
    # Notice how we are deliberately skipping the first and last items in the list
    # because we don't need to compare them against themselves
    for i in range(1, len(fp_list)):
        # Compare the current fingerprint against all the previous ones in the list
        similarities = DataStructs.BulkTanimotoSimilarity(fp_list[i], fp_list[:i])
        # Since we need a distance matrix, calculate 1-x for every element in similarity matrix
        dissimilarity_matrix.extend([1 - x for x in similarities])
    return dissimilarity_matrix


def cluster_fingerprints(fingerprints, cutoff):
    """Cluster fingerprints
    Parameters:
        fingerprints
        cutoff: threshold for the clustering
    """
    # Calculate Tanimoto distance matrix
    distance_matrix = tanimoto_distance_matrix(fingerprints)
    # Now cluster the data with the implemented Butina algorithm:
    clusters = Butina.ClusterData(distance_matrix, len(fingerprints), cutoff, isDistData=True)
    clusters = sorted(clusters, key=len, reverse=True)
    return clusters



# function to get properties from a list of CIDs

def get_properties_from_CIDs(CIDs):
    """Get molecular properties for a list of CIDs from PubChem. Returns DataFrame or 'error' on failure."""

    # batch processing of CIDs
    chunk_size = 400
    all_data = []

    for i in range(0, len(CIDs), chunk_size):
        chunk = CIDs[i:i+chunk_size]
        cid_string = ','.join(map(str, chunk))
        
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_string}/property/Title,SMILES,XLogP,HBondDonorCount,TPSA,HBondAcceptorCount,MolecularWeight,MolecularFormula/csv'

        try:
            response = requests.get(url, timeout=10)
        except requests.RequestException as e:
            return "error"

        if response.status_code != 200:
            return "error"

        try:
            data = pd.read_csv(StringIO(response.text))
        except Exception as e:
            return "error"

        if data.empty:
            return "error"
        
        all_data.append(data)
        time.sleep(0.3)  # Be polite to PubChem
    
    return pd.concat(all_data, ignore_index=True)



def get_synonyms_from_CID(CID):
    """Get up to 10 synonyms from a PubChem CID. Returns 'error' if the request fails or no synonyms are found."""

    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{CID}/synonyms/TXT'

    try:
        response = requests.get(url, timeout=10)
    except requests.RequestException as e:
        return "error"

    if response.status_code != 200:
        return "error"

    synonyms = response.text.strip().split('\n')[:10]

    if not synonyms or synonyms == ['']:
        return "error"

    # Return up to 10 joined synonyms
    joined_synonym = ";".join(synonyms)
    return joined_synonym



# function to get CID from smiles

def get_CID_from_SMILES(smiles_b64):
    """Get PubChem CID from base64-encoded SMILES. Returns CID string or 'error'."""

    try:
        smiles = base64.b64decode(smiles_b64).decode('utf-8').strip()
    except Exception:
        return "error"

    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/TXT"
    headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    data = {'smiles': smiles}

    try:
        response = requests.post(url, data=data, headers=headers, timeout=10)
    except requests.RequestException:
        return "error"

    if response.status_code != 200:
        return "error"

    cid = response.text.strip()

    if not cid or cid == "0":
        return "error"

    return cid


def get_cid_from_dtxsid(dtxsid):
    if not dtxsid or not isinstance(dtxsid, str):
        return "error"

    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/xref/RegistryID/{dtxsid}/cids/TXT'

    try:
        response = requests.get(url, timeout=10)
    except requests.RequestException:
        return "error"

    if response.status_code != 200:
        return "error"

    cid = response.text.strip()

    if not cid or not cid.isdigit():
        return "error"

    return cid


def check_cid(CID):
    

    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{CID}/property/Title/TXT'

    try:
        response = requests.get(url, timeout=10)
    except requests.RequestException:
        return "error"

    if response.status_code != 200:
        return "error"

    name = response.text.strip()

    if not name:
        return "error"

    return name


### function to integrate silhouette score plots to app


def create_full_distance_matrix(fp_list):
    """this create the full matrix; upper and lower triangle for silhoutte score calculations
     sklearn needs the full matrix which we do not create for Butina where its a flat array

    """
    n = len(fp_list)
    distance_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i+1, n):
            # Calculate Tanimoto similarity
            similarity = DataStructs.TanimotoSimilarity(fp_list[i], fp_list[j])
            distance = 1 - similarity
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance  # Make symmetric
    
    return distance_matrix



def clusters_to_labels(clusters, n_molecules):
    """Convert Butina cluster format to label array for silhouette analysis"""
    labels = np.full(n_molecules, -1)  # Initialize with -1 (noise/unclustered)
    
    for cluster_id, cluster in enumerate(clusters):
        for molecule_idx in cluster:
            labels[molecule_idx] = cluster_id
    
    return labels


def calculate_silhouette_scores(fingerprints, cutoff_range=None, step=0.05):
    """
    Calculate silhouette scores for different distance cutoffs
    
    Parameters:
        fingerprints: list of molecular fingerprints
        cutoff_range: tuple (min_cutoff, max_cutoff) or None for automatic range
        step: step size for cutoff values
    
    Returns:
        dict with cutoffs as keys and silhouette scores as values
    """
    if cutoff_range is None:
        cutoff_range = (0.1, 0.9)  # Default range
    
    cutoffs = np.arange(cutoff_range[0], cutoff_range[1] + step, step)
    silhouette_scores = {}
    cluster_counts = {}
    valid_cutoffs = []
    
    # Create full distance matrix once (expensive operation)
    print("Creating full distance matrix...")
    full_distance_matrix = create_full_distance_matrix(fingerprints)
    
    print(f"Testing {len(cutoffs)} different cutoff values...")
    
    for cutoff in cutoffs:
        try:
            # Get clusters for this cutoff
            clusters = cluster_fingerprints(fingerprints, cutoff)
            
            # Convert to labels format
            labels = clusters_to_labels(clusters, len(fingerprints))
            
            # Skip if we have less than 2 clusters or too many singletons
            unique_labels = np.unique(labels)
            n_clusters = len(unique_labels[unique_labels >= 0])  # Exclude -1 (noise)
            
            if n_clusters < 2:
                print(f"Cutoff {cutoff:.3f}: Only {n_clusters} cluster(s), skipping...")
                continue
            
            # Calculate silhouette score
            sil_score = silhouette_score(full_distance_matrix, labels, metric='precomputed')
            silhouette_scores[cutoff] = sil_score
            cluster_counts[cutoff] = n_clusters
            valid_cutoffs.append(cutoff)
            
            print(f"Cutoff {cutoff:.3f}: {n_clusters} clusters, silhouette = {sil_score:.3f}")
            
        except Exception as e:
            print(f"Error at cutoff {cutoff:.3f}: {e}")
            continue
    
    return silhouette_scores, cluster_counts, valid_cutoffs



# Main execution function
def run_silhouette_analysis(fingerprints, cutoff_range=(0.1, 0.9), step=0.05):
    """
    Complete silhouette analysis workflow
    
    Parameters:
        fingerprints: list of molecular fingerprints
        cutoff_range: tuple (min_cutoff, max_cutoff)
        step: step size for cutoff values
        save_plot: path to save plot (optional)
    
    Returns:
        tuple: (best_cutoff, silhouette_scores, cluster_counts)
    """
    print(f"Starting silhouette analysis for {len(fingerprints)} molecules...")
    
    # Calculate silhouette scores
    silhouette_scores, cluster_counts, valid_cutoffs = calculate_silhouette_scores(
        fingerprints, cutoff_range, step
    )
    
    if not silhouette_scores:
        print("No valid clustering results found. Try adjusting the cutoff range.")
        return None, None, None
    
    
    return silhouette_scores


# Main execution function
def get_cluster_counts(fingerprints, cutoff_range=(0.1, 0.9), step=0.05):
    """
    Complete silhouette analysis workflow
    
    Parameters:
        fingerprints: list of molecular fingerprints
        cutoff_range: tuple (min_cutoff, max_cutoff)
        step: step size for cutoff values
        save_plot: path to save plot (optional)
    
    Returns:
        tuple: (best_cutoff, silhouette_scores, cluster_counts)
    """
    print(f"Starting silhouette analysis for {len(fingerprints)} molecules...")
    
    # Calculate silhouette scores
    silhouette_scores, cluster_counts, valid_cutoffs = calculate_silhouette_scores(
        fingerprints, cutoff_range, step
    )
    
    if not silhouette_scores:
        print("No valid clustering results found. Try adjusting the cutoff range.")
        return None, None, None
    
    
    return cluster_counts








##################################################################################
## extra code beyond this point



# Function to get CID from CASRN
# do not use this function in app but might need it later
def get_CID_from_CAS(CASRN):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/xref/RegistryID/{CASRN}/cids/TXT'
    response = requests.get(url)
    return response.text.strip()


