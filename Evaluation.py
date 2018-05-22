from __future__ import division
import numpy as np
from chemfp import search

DEBUG = False

# calculate number of actives and decoys in a cluster
# for numpy molecular representations
def getNumberOfMolecules(cluster, activeSubSet, decoySubSet):
    # get the number of molecules in cluster which belong to the subset given
    nrows, ncols = activeSubSet.shape
    dtype = {'names': ['f{}'.format(i) for i in range(ncols)],
             'formats': ncols * [activeSubSet.dtype]}
    # print dtype
    activesInCluster = np.intersect1d(activeSubSet.view(dtype), cluster.view(dtype))
    decoysInCluster = np.intersect1d(decoySubSet.view(dtype), cluster.view(dtype))
    return len(activesInCluster), len(decoysInCluster)


# return number of actives and Decoys in cluster - with molecules
# being represented as Fingerprints
def getNumberOfMoleculesFP(cluster, activeSubSet, comparison_type):
    actives_in_cluster = 0
    totalNumInCluster = len(cluster)
    for clusterEl in cluster:
        if comparison_type == "RDKIT":
            if clusterEl in activeSubSet:
                #print clusterEl
                actives_in_cluster = actives_in_cluster + 1
        elif comparison_type == "CHEMFP_ID":
            active_mols = search.contains_fp(clusterEl[1], activeSubSet)
            #active_mol = activeSubSet.get_by_id(id=clusterEl)
            #print len(active_mols)
            if len(active_mols) > 0:
                #print active_mol[0]
                #print "Found active molecule"
                actives_in_cluster = actives_in_cluster + 1

    if DEBUG:
        print len(cluster)
        print "Actives in cluster = ", actives_in_cluster, " decoys ", totalNumInCluster - actives_in_cluster
    return actives_in_cluster, (totalNumInCluster - actives_in_cluster)


def calculatePercentage(num1, num2):
    return (num1 / (num1 + num2)) * 100


# calculation of QPI having labels(label represents clusterid)
def calculateQPIWithLabels(labels, num_clusters, mergedFingerprints, activeDescriptorsSubset,comparison_type):
    clusters = []
    for x in range(0, num_clusters):
        indeces = np.where(labels == (x))[0]
        cluster = [mergedFingerprints[i] for i in indeces]
        clusters.append(cluster)
    return calculate_qpi(clusters, activeDescriptorsSubset, len(mergedFingerprints),comparison_type)

# calculate QPI passing as parameter, clusters with index of each molecule in cluster
def calculateQPIWithIndeces(indeces_clusters, mergedFingerprints, activeDescriptorsSubset,comparison_type):
    clusters = []
    for indeces in indeces_clusters:
        if comparison_type == "RDKIT":
            cluster = [mergedFingerprints[i] for i in indeces]
            clusters.append(cluster)
        elif comparison_type == "CHEMFP_ID":
            #merged fingerprints is an arena in chemfp
            #print [mergedFingerprints[i] for i in indeces]
            cluster = [mergedFingerprints[i] for i in indeces]
            clusters.append(cluster)
    return calculate_qpi(clusters, activeDescriptorsSubset, len(mergedFingerprints), comparison_type)


# calculation of QPI passing the clusters as parameters
def calculate_qpi(clusters, active_descriptors, total_number_molecules, comparison_type):
    percentage_actives_in_set = calculatePercentage(len(active_descriptors), total_number_molecules)
    p = 0
    q = 0
    r = 0
    s = 0
    total_actives = 0
    for cluster in clusters:
        # numActives, numDecoys = getNumberOfMolecules(cluster, activeDescriptorsSubset, decoyDescriptorsSubset)
        num_actives, num_decoys = getNumberOfMoleculesFP(cluster, active_descriptors, comparison_type)
        #print "Actives in cluster: ", num_actives, " Decoys in cluster: ", num_decoys
        if DEBUG:
            print "Actives ", num_actives, " Decoys ", num_decoys
        total_actives = total_actives + num_actives
        percentage_actives_cluster = calculatePercentage(num_actives, num_decoys)

        if (num_actives == 1) and (num_decoys == 0):
            s += 1
        elif percentage_actives_cluster > percentage_actives_in_set:
            p += num_actives
            q += num_decoys
        else:
            r += num_actives
    qpi = p / (p + q + r + s)

    print "p", p , " q " , q,  " r ", r , " s ", s
    print "Total Actives " ,total_actives
    return qpi

#region PRECISION RECALL F1

def calculate_f1(clusters, active_fps):
    a = 0
    n = 0
    A = len(active_fps)
    total_actives = 0
    precision = 0
    recall = 0
    fscore = 0
    results = []

    max_fscore = 0
    max_result = []
    cluster_id = 0
    for cluster in clusters:
        n = len(cluster)
        num_actives, num_decoys = getNumberOfMoleculesFP(cluster, active_fps, "RDKIT")
        total_actives += num_actives
        a = num_actives

        if DEBUG:
            print "Actives ", num_actives, " Decoys ", num_decoys
        precision = a/n
        recall = a/A

        if precision != 0 and recall != 0:
            fscore = (2*precision * recall)/(precision + recall)
        else:
            fscore = 0

        if fscore > max_fscore:
            max_fscore = fscore
            max_result = [precision, recall, fscore]

        results.append([cluster_id, precision, recall, fscore, num_actives, num_decoys])
        cluster_id += 1

    print "Total Actives ", total_actives
    return max_result, results
#endregion PRECISION RECALL F1

#region ENRICHMENT FACTOR
def calculate_enrichment(clusters, active_fps):
    dataset_actives = len(active_fps)
    dataset_size = 0

    for cluster in clusters:
        dataset_size += len(cluster)

    print "dataset actives, ", dataset_actives
    print "dataset size, ", dataset_size
    results = []

    max_ef = 0
    cluster_id = 0
    for cluster in clusters:
        cluster_size = len(cluster)
        num_actives, num_decoys = getNumberOfMoleculesFP(cluster, active_fps, "RDKIT")

        if DEBUG:
            print "Actives ", num_actives, " Decoys ", num_decoys

        enrichment = (num_actives/dataset_actives)/(cluster_size/dataset_size)

        if enrichment > max_ef:
            max_ef = enrichment

        results.append([cluster_id, enrichment, num_actives, num_decoys])
        cluster_id += 1

    return max_ef, results
#endregion ENRICHMENT FACTOR