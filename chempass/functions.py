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
    """Get molecular properties for a list of CIDs from PubChem. Returns DataFrame or 'error' on failure."""

    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{CIDs}/property/Title,SMILES,XLogP,HBondDonorCount,TPSA,HBondAcceptorCount,MolecularWeight,MolecularFormula/csv'

    try:
        response = requests.get(url, timeout=10)
    except requests.RequestException as e:
        #print(f"❌ Request error for CIDs {CIDs}: {e}")
        return "error"

    if response.status_code != 200:
        #print(f"⚠️ PubChem returned HTTP {response.status_code} for CIDs {CIDs}")
        return "error"

    try:
        data = pd.read_csv(StringIO(response.text))
    except Exception as e:
        #print(f"❌ Failed to parse CSV for CIDs {CIDs}: {e}")
        return "error"

    if data.empty:
        #print(f"⚠️ No property data returned for CIDs {CIDs}")
        return "error"

    return data



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





##################################################################################
## extra code beyond this point



# Function to get CID from CASRN
# do not use this function in app but might need it later
def get_CID_from_CAS(CASRN):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/xref/RegistryID/{CASRN}/cids/TXT'
    response = requests.get(url)
    return response.text.strip()


