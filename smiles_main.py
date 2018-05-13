from rdkit import Chem
from MolecularRepresentations import *
from Clustering import *
from Evaluation import *
from Reader import *
import os


BUTINA = True
KMEANS = False
WARDS = False
RDKIT = False
JPCLUSTERING = False
LEADER = True
CONVERT = False

def sfd_to_smiles(input_file, output_file):
    mols = [mol for mol in Chem.SDMolSupplier(input_file) if mol != None]

    #print os.getcwd()
    #os.chdir("SmilesActives")
    writer = SmilesWriter(output_file)
    for mol in mols:
        writer.write(mol)
    writer.close()

if __name__ == "__main__":
    file = '../mols/merged/1000ABL1.smi'

    if BUTINA:
        fingerprints, molecules = read_smiles_mol_fps(file)

        tuple_list, neighbours_list = get_neighbours_list(fingerprints, 0.3)
        clusters = ButinaClustering(tuple_list, neighbours_list, len(fingerprints))

        mol_clusters = change_indeces_to_smiles(clusters,molecules)
        output_cluster_results(mol_clusters)
        #for cluster in clusters:
        #    print cluster

    #Convert sdf files to SMILES
    if CONVERT:
        actives_list = ['ABL1', 'AMPC', 'ANDR', 'FNTA', 'Renin', 'THB', 'THRB']
        for active_file in actives_list:
            sfd_to_smiles('dataset/'+active_file+'/actives_final.sdf','dataset/SmilesActives/'+active_file+'.smi')