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
import matplotlib.pyplot as plt
from rdkit.Chem import Draw
import requests
import re
from io import StringIO
import base64



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
    """Function to take in valid CIDs (no missing CIDs are expected to be present in this input list) in batch and provide molecular properties"""
    
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{CIDs}/property/Title,SMILES,XLogP,HBondDonorCount,TPSA,HBondAcceptorCount,MolecularWeight,MolecularFormula/csv'
    response = requests.get(url)

    data = pd.read_csv(StringIO(response.text))
    
    return data




# Function to get CID from CASRN
# do not use this function in app but might need it later
def get_CID_from_CAS(CASRN):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/xref/RegistryID/{CASRN}/cids/TXT'
    response = requests.get(url)
    return response.text.strip()

# function to get CID from DTXSID
def get_cid_from_dtxsid(dtxsid):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/xref/RegistryID/{dtxsid}/cids/TXT'
    response = requests.get(url)
    if  response.status_code == 200:
        if response.text.strip():
            return response.text.strip()
        else:
            return 'error'          
    else:
        return 'error'


def get_synonyms_from_CID(CID):
    """Function to take in valid CIDs (no missing CIDs are expected to be present in this input list) one by one and provide upto ten synonyms"""


    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{CID}/synonyms/TXT'
    response = requests.get(url)

    synonyms = response.text.strip().split('\n')[:10]
    joined_synonym = ";".join(synonyms)

    return joined_synonym



# function to get CID from smiles
def get_CID_from_SMILES(smiles_b64):
    """Function to take in valid SMILES one by one and provide CIDs. This function is set to handle special characters in R like back slashes or dots"""


    smiles = base64.b64decode(smiles_b64).decode('utf-8')
    
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/TXT"

    headers = {
        'Content-Type': 'application/x-www-form-urlencoded'
    }

    data = {
        'smiles': smiles
    }

    response = requests.post(url, data=data, headers=headers)

    if  response.status_code == 200:
        if response.text.strip():
            return response.text.strip()
        else:
            return 'error'          
    else:
        return 'error'

   





