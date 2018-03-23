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
RDKIT = False


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

if __name__ == "__main__":
    #RdKitButina()

    percentage_list = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
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
            activeDescriptors = read_fingerprints('dataset/actives_final.sdf')
            decoyDescriptors = read_fingerprints('dataset/decoys_final.sdf')
            for perc in percentage_list:
                activeDescriptorsSubset, decoyDescriptorsSubset = getPercentageMolecules(perc, activeDescriptors, decoyDescriptors)

                print "Clustering ", len(activeDescriptorsSubset), " active molecules and ", len(decoyDescriptorsSubset), " decoy molecules"
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

    #        if not RDKIT:
    #            arena_actives = read_chemfp("dataset/actives_final.sdf")
    #            arena_all = read_chemfp("dataset/merged.sdf")


    #            arena_actives_subset, arena_subset = get_arena_percentage(0.8, arena_actives, arena_all)
    #            print "Active Subset Size ", len(arena_actives_subset)
    #            print "Total Subset Size ", len(arena_subset)

    #            cluster_result = taylor_butina_cluster(arena_subset, 0.8)
    #            clusters = combine_clusters(cluster_result)
    #            qpi = calculate_qpi(clusters, arena_actives, len(arena_subset), "CHEMFP_ID")
    #            print "Qpi for clustering Butina: ", qpi

            #count = 0
            #print len(arena_all)
            #print len(arena_actives)
            #for mol in arena_all:
            #    mol_id = mol[0]
            #    # print mol_id
            #    active_mol = arena_actives.get_by_id(id=mol_id)
            #    if active_mol is not None:
            #        print "Found active ", mol
            #        count = count + 1
            #print "found ", count, " actives "
    if KMEANS:
        # KMEANS Clustering
        labels = KMeansClustering(mergedFingerprints)
        qpi = calculateQPIWithLabels(labels, 10, mergedFingerprints, activeDescriptorsSubset, decoyDescriptorsSubset)
        print "Qpi for clustering KMEANS: ", qpi

    if WARDS:
        # Wards Clustering
        c_tree = WardsClustering(dists, len(mergedFingerprints))

        ward_clusters = GetHierarchicalLevel(c_tree[0], 3)
        qpi = calculateQPIWithIndeces(ward_clusters, mergedFingerprints, activeDescriptorsSubset, decoyDescriptorsSubset)
        print "Qpi for clustering WARD level 3 : ", qpi

        ward_clusters = GetHierarchicalLevel(c_tree[0], 4)
        qpi = calculateQPIWithIndeces(ward_clusters, mergedFingerprints, activeDescriptorsSubset, decoyDescriptorsSubset)
        print "Qpi for clustering WARD level 4 : ", qpi

        ward_clusters = GetHierarchicalLevel(c_tree[0], 5)
        qpi = calculateQPIWithIndeces(ward_clusters, mergedFingerprints, activeDescriptorsSubset, decoyDescriptorsSubset)
        print "Qpi for clustering WARD level 5 : ", qpi

