'''
This file contains method related to file reading
'''
import chemfp
from rdkit import Chem

from MolecularRepresentations import *

def read_chemfp(input_file):
    reader = chemfp.read_molecule_fingerprints("RDKit-Morgan fpSize=1024", input_file)
    arena_all = chemfp.load_fingerprints(reader, reorder=False)
    return arena_all

def read_fingerprints(input_file):
    molecules = Chem.SDMolSupplier(input_file)

    moleculeFingerprints = createFingerprint(molecules)
    return moleculeFingerprints

def read_smiles_fingerprints(input_file):
    molecules = Chem.SmilesMolSupplier(input_file)

    moleculeFingerprints = createFingerprint(molecules)
    return moleculeFingerprints

def read_smiles_mol_fps(input_file):
    molecules = Chem.SmilesMolSupplier(input_file)

    moleculeFingerprints = createFingerprint(molecules)
    return moleculeFingerprints, molecules

def read_smiles_fps_clusters(input_file):
    molecules = Chem.SmilesMolSupplier(input_file)
    #list_Res = molecules[0].GetPropNames()

    clusters = create_cluster_fingerprints(molecules)
    return clusters
    #for mol in molecules:
    #   print mol.GetProp('Cluster')