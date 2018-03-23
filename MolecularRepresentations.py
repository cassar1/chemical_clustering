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