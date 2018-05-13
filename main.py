from rdkit import Chem
from MolecularRepresentations import *
from Clustering import *
from Evaluation import *
from Reader import *
import chemfp
from chemfp import search
from itertools import groupby

BUTINA = False
KMEANS = False
WARDS = False
RDKIT = True
JPCLUSTERING = False
LEADER = True


def RdKitButina():
    activeDescriptors = read_fingerprints('dataset/actives_final.sdf')
    decoyDescriptors = read_fingerprints('dataset/decoys_final.sdf')
    activeDescriptorsSubset, decoyDescriptorsSubset = getPercentageMolecules(0.2, activeDescriptors, decoyDescriptors)

    mergedFingerprints = activeDescriptorsSubset[:]
    mergedFingerprints.extend(decoyDescriptorsSubset)
    dists, dists_2d = CalcSimilarities(mergedFingerprints)
    clusters = ButinaClusteringOriginal(dists, len(mergedFingerprints))
    # for key, group in groupby(clusters, key=lambda x: len(x)):
    #    print key, len(list(group))
    qpi = calculateQPIWithIndeces(clusters, mergedFingerprints, activeDescriptorsSubset, "RDKIT")
    print "Qpi for Butina: ", qpi
    print "-------------------------------------------------------------------------------------------"


def jp_clustering_setup(percentage_mol):
    #https://iwatobipen.wordpress.com/2016/03/13/jp-clustering-with-python-and-rdkit/
    activeDescriptors = read_fingerprints('dataset/THB/actives_final.sdf')
    decoyDescriptors = read_fingerprints('dataset/AMPC/actives_final.sdf')
    activeDescriptorsSubset, decoyDescriptorsSubset = getPercentageMolecules(percentage_mol, activeDescriptors, decoyDescriptors)
    print "Clustering ", len(activeDescriptorsSubset), " active molecules and ", len(decoyDescriptorsSubset), " decoy molecules"
    mergedFingerprints = activeDescriptorsSubset[:]
    mergedFingerprints.extend(decoyDescriptorsSubset)
    clusters = jp_clustering(mergedFingerprints, 16, 4)

    qpi = calculate_qpi(clusters, activeDescriptorsSubset, len(mergedFingerprints),  "RDKIT")
    print "Qpi for JP Clustering: ", qpi
    print "-------------------------------------------------------------------------------------------"

def leader_clustering_setup(actives_file, decoys_file, percentage_mol, threshold):
    activeDescriptors = read_fingerprints(actives_file)
    decoyDescriptors = read_fingerprints(decoys_file)
    activeDescriptorsSubset, decoyDescriptorsSubset = getPercentageMolecules(percentage_mol, activeDescriptors, decoyDescriptors)
    print "Clustering ", len(activeDescriptorsSubset), " active molecules and ", len(decoyDescriptorsSubset), " decoy molecules"
    mergedFingerprints = activeDescriptorsSubset[:]
    mergedFingerprints.extend(decoyDescriptorsSubset)
    clusters = leader_algorithm(mergedFingerprints, threshold)

    qpi = calculate_qpi(clusters, activeDescriptorsSubset, len(mergedFingerprints),  "RDKIT")
    print "Qpi for Leader Clustering: ", qpi
    print "-------------------------------------------------------------------------------------------"


if __name__ == "__main__":
    #RdKitButina()
    actives_file = 'dataset/THB/actives_final.sdf'
    decoys_file = 'dataset/THB/decoys_final.sdf'
    percentage_list = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
    #percentage_list = [1]
    if BUTINA:
        if not RDKIT:
            print "CHEMFP"
            arena_actives = read_chemfp("dataset/actives_final.sdf")
            arena_all = read_chemfp("dataset/merged.sdf")

            for perc in percentage_list:
                print "Clustering for ", perc, " molecules "
                arena_actives_subset, arena_subset = get_arena_percentage(perc, arena_actives, arena_all)
                start_time = time.time()
                print "Actives Subset length ", len(arena_actives_subset)
                print "Merged Subset length " ,len(arena_subset)
                similarity_table = search.threshold_tanimoto_search_symmetric(
                    arena_subset, threshold=0.5)

                #centroid_table = sorted(((len(indices), i, indices)
                #                         for (i, indices) in enumerate(similarity_table.iter_indices())),
                #                        reverse=True)
                tuple_list = sorted(((len(indices), i)
                                 for (i, indices) in enumerate(similarity_table.iter_indices())),
                                reverse=True)
                neighbours_list = [indices for (i, indices) in enumerate(similarity_table.iter_indices())]
                print "time taken to calculate neighbours: ", time.time() - start_time
                clusters = ButinaClustering(tuple_list, neighbours_list, len(arena_subset))

                #for key,group in groupby(clusters, key=lambda x:len(x)):
                #    print key, len(list(group))

                #qpi = calculate_qpi(clusters, arena_actives, len(arena_subset), "CHEMFP_ID")
                qpi = calculateQPIWithIndeces(clusters, arena_subset, arena_actives_subset, "CHEMFP_ID")
                print "Qpi for clustering Butina: ", qpi
                print "-----------------------------------------------------------------------------------"
            #print "RDKIT"
        if RDKIT:
            activeDescriptors = read_fingerprints(actives_file)
            decoyDescriptors = read_fingerprints(decoys_file)
            for perc in percentage_list:
                activeDescriptorsSubset, decoyDescriptorsSubset = getPercentageMolecules(perc, activeDescriptors, decoyDescriptors,shuffle=True)

                print "Clustering ",perc ," percentage ", len(activeDescriptorsSubset), " active molecules and ", len(decoyDescriptorsSubset), " decoy molecules"
                mergedFingerprints = activeDescriptorsSubset[:]
                mergedFingerprints.extend(decoyDescriptorsSubset)
                print "Calculating neighbours"
                #dists, dists_2d = CalcSimilarities(mergedFingerprints)
                tuple_list, neighbours_list = get_neighbours_list(mergedFingerprints, 0.8)
                print "Finished calculating neighbours"

                # Butina Clustering
                clusters = ButinaClustering(tuple_list, neighbours_list, len(mergedFingerprints))
                #for key, group in groupby(clusters, key=lambda x: len(x)):
                #    print key, len(list(group))
                qpi = calculateQPIWithIndeces(clusters,mergedFingerprints, activeDescriptorsSubset, "RDKIT")
                print "Qpi for Butina: ", qpi
                print "-----------------------------------------------------------------------------------"

    if KMEANS:
        # KMEANS Clustering
        labels = KMeansClustering(mergedFingerprints)
        qpi = calculateQPIWithLabels(labels, 10, mergedFingerprints, activeDescriptorsSubset, decoyDescriptorsSubset)
        print "Qpi for clustering KMEANS: ", qpi

    if WARDS:
        #percentage_mols = 0.01
        if not RDKIT:
            print "CHEMFP"
            arena_actives = read_chemfp(actives_file)
            arena_all = read_chemfp(decoys_file)
            arena_actives_subset, arena_subset = get_arena_percentage(0.1, arena_actives, arena_all)

            dists = distance_matrix_1d(arena_subset)
            #dists2 = distance_matrix(arena_subset)
            c_tree = WardsClustering(dists, len(arena_subset))

            for i in range(1,len(arena_subset),100):
                ward_clusters = get_hierarchical_level(c_tree[0], i)
                qpi = calculateQPIWithIndeces(ward_clusters, arena_subset, arena_actives_subset, "CHEMFP_ID")
                print "Qpi for clustering WARD level ", i ," : ", qpi
            #c_tree = WardsClusteringScipy(dists, len(arena_subset))
            #OutputDendrogram(c_tree)
        if RDKIT:
            SKLEARN = False
            for perc in percentage_list:
                activeDescriptors = read_fingerprints(actives_file)
                decoyDescriptors = read_fingerprints(decoys_file)
                activeDescriptorsSubset, decoyDescriptorsSubset = getPercentageMolecules(perc, activeDescriptors,
                                                                                         decoyDescriptors)

                print "Clustering ", len(activeDescriptorsSubset), " active molecules and ", len(
                    decoyDescriptorsSubset), " decoy molecules"
                mergedFingerprints = activeDescriptorsSubset[:]
                mergedFingerprints.extend(decoyDescriptorsSubset)

                if SKLEARN:
                    num_clusters = len(mergedFingerprints)
                    for clusters in range(50,len(mergedFingerprints),10):
                        print clusters, " clusters"
                        fingerprintsToNPArr(mergedFingerprints)
                        labels = ward_clustering_sklearn(clusters, mergedFingerprints)

                        qpi = calculateQPIWithLabels(labels, clusters, mergedFingerprints, activeDescriptorsSubset, "RDKIT")
                        print "Qpi for clustering Wards: ", qpi
                else:
                    dists = calc_distance_1d(mergedFingerprints)

                    c_tree = WardsClustering(dists, len(mergedFingerprints))
                    print "Finished clustering"

                    # Wards Clustering
                    if perc == 2:
                        i = 10

                        #for i in range(2,len(mergedFingerprints),10):
                        while i < len(mergedFingerprints):
                            print i
                            ward_clusters = get_hierarchical_level1(c_tree[0], i)
                            qpi = calculateQPIWithIndeces(ward_clusters, mergedFingerprints, activeDescriptorsSubset, "RDKIT")
                            if qpi > 0.5:
                                print "Qpi for clustering WARD level ", i ," with ", len(ward_clusters), " clusters : ", qpi
                                print "-----------------------------------------------------------------------------------------"
                                i += 10
                            else:
                                i += 10
    if JPCLUSTERING:
        for perc in percentage_list:
            jp_clustering_setup(perc)

    if LEADER:
        for perc in percentage_list:
            leader_clustering_setup(actives_file, decoys_file, perc, 0.5)


