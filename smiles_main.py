from rdkit import Chem
from MolecularRepresentations import *
from Clustering import *
from Evaluation import *
from Reader import *
import os


BUTINA = False
WARDS = True
JPCLUSTERING = False
LEADER = False
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
    #file = '../mols/compounds5.smi'
    #file = '../mols/merged/1000ABL1.smi'
    file = 'dataset/SmilesMerged/ABL1merged.smi'

    fingerprints, molecules = read_smiles_mol_fps(file)
    similarity_list = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

    if BUTINA:
        for similarity_threshold in similarity_list:
            print "Calculating Butina at similarity ", similarity_threshold
            tuple_list, neighbours_list = get_neighbours_list(fingerprints, similarity_threshold)
            #for tpl in tuple_list:
            #    print tpl[1],",", neighbours_list[tpl[1]]
            clusters = ButinaClustering(tuple_list, neighbours_list, len(fingerprints))
            clusters.sort(key=len)
            #for cluster in clusters:
            #    print cluster
            mol_clusters = change_indeces_to_smiles(clusters,molecules)
            output_cluster_results(mol_clusters, 'butinaRes' + str(similarity_threshold), 'butina')
    if WARDS:
        dists = calc_distance_1d(fingerprints)

        c_tree = WardsClustering(dists, len(fingerprints))
        print "Finished clustering"

        # Wards Clustering
        i = 4000
        while i < len(fingerprints):
            ward_clusters = get_hierarchical_level(c_tree[0], i)
            mol_clusters = change_indeces_to_smiles(ward_clusters, molecules)
            output_cluster_results(mol_clusters,'wardResult' + str(i), 'wards')
            i += 100

    if JPCLUSTERING:
        #www.daylight.com/dayhtml/doc/cluster/
        #suggests 16 as default size
        neighbours = 16
        k_min_neighbours = 4
        #,(16,8),(16,12)
        parameters = [(16,4)]
        #clusters = jp_clustering(fingerprints, neighbours, k_min_neighbours)
        cluster_gen = jp_clustering_neighbours(molecules)

        for param in parameters:
            print "Grouping ",param[0], " ", param[1]
            mol_clusters = jp_clustering_results(cluster_gen, param[0], param[1])

            #mol_clusters = change_indeces_to_smiles(clusters, molecules)
            output_cluster_results(mol_clusters, 'jarvisResult_'+str(param[0]) + '_' + str(param[1]), 'jarvis_pat')

    if LEADER:
        for similarity_threshold in similarity_list:
            clusters = leader_algorithm_ind(fingerprints, similarity_threshold)
            mol_clusters = change_indeces_to_smiles(clusters, molecules)
            output_cluster_results(mol_clusters, 'leader' + str(similarity_threshold), 'leaders')

    #Convert sdf files to SMILES
    if CONVERT:
        convert_actives = False
        convert_merged = True
        if convert_actives:
            actives_list = ['ABL1', 'AMPC', 'ANDR', 'FNTA', 'Renin', 'THB', 'THRB']
            for active_file in actives_list:
                sfd_to_smiles('dataset/'+active_file+'/actives_final.sdf','dataset/SmilesActives/'+active_file+'.smi')
        if convert_merged:
            merged_list = ['ABL1', 'AMPC', 'ANDR', 'FNTA', 'Renin', 'THB', 'THRB']
            for merged_file in merged_list:
                sfd_to_smiles('dataset/' + merged_file + '/merged.sdf',
                              'dataset/SmilesMerged/' + merged_file + 'merged.smi')

