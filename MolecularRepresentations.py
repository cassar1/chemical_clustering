from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit import Chem

from rdkit import DataStructs
from rdkit.ML.Descriptors import MoleculeDescriptors
import numpy as np

# gets fingerprint - which contains counts
# AllChem.GetMorganFingerprint
# gets fingerprint - in bit vector format (true/false)
# AllChem.GetMorganFingerprintAsBitVect
def createFingerprint(mols, counts= False):
    if not counts:
        fingerprints = [AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024) for m in mols]
        return fingerprints
    else:
        info = {}
        fingerprints = [AllChem.GetMorganFingerprint(m,2,bitInfo=info) for m in mols]
        return fingerprints

def readAndCreateFingerprint(file_name, counts=False):
    if not counts:
        fingerprints = [AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024) for m in
                        Chem.ForwardSDMolSupplier(file_name, removeHs=False) if m is not None]
        return fingerprints
    else:
        info = {}
        fingerprints = [AllChem.GetMorganFingerprint(m, 2, bitInfo=info) for m in
                        Chem.ForwardSDMOLSupplier(file_name, removeHs=False) if m is not None]
        return fingerprints

# MACCS Fingerprints
def createMACCSFingerprint(mols):
    fingerprints = [MACCSkeys.GenMACCSKeys(m) for m in mols]
    for fingerprint in fingerprints:
        print fingerprint.ToBitString()
    return fingerprints


# Convert bit vector fingerprint to np array
# does not work with count based fingerprints
def fingerprintsToNPArr(fps):
    # print fps
    np_fps = []
    for fp in fps:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        np_fps.append(arr)
    return np_fps


def create_cluster_fingerprints(mols, counts= False):
    clusters = {}

    if not counts:
        for mol in mols:
            cluster_id = mol.GetProp('Cluster')
            if cluster_id not in clusters:
                clusters[cluster_id] = []
            clusters[cluster_id].append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024))
    else:
        for mol in mols:
            cluster_id = mol.GetProp('Cluster')
            if cluster_id not in clusters:
                clusters[cluster_id] = []
            clusters[cluster_id].append(AllChem.GetMorganFingerprint(mol,2))

    return clusters.values()