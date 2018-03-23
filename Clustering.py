import time
from numpy.f2py.auxfuncs import throw_error
from rdkit.ML.Cluster import Butina, Murtagh
from sklearn.cluster import KMeans
from MolecularRepresentations import *
from sklearn.cluster import AgglomerativeClustering

from Helpers import *
import chemfp
from chemfp import search

# Butina Clustering Algorithm
def ButinaClusteringOriginal(dists, nfps):
    print "-------------------------------------------------"
    print "starting Butina clustering"
    # now cluster the data:
    start_time = time.time()
    cs = Butina.ClusterData(dists,nfps,0.7,isDistData=True, reordering=True)

    print "time taken: ", time.time() - start_time
    return cs

def ButinaClustering(tuple_list, neighbours_list, nfps):
    print "-------------------------------------------------"
    print "starting Butina clustering"
    # now cluster the data:
    start_time = time.time()
    #cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True, reordering=True)
    cs = custom_butina(tuple_list, neighbours_list, nfps, reordering=False)
    print "time taken: ", time.time() - start_time
    return cs


# Butina Clustering without calculation of Neighbours
# copied from RDKIT Library
def custom_butina(tLists , nbrLists, nPts, reordering=False):
    res = []
    seen = [0] * nPts
    while tLists:
        _, idx = tLists.pop(0)
        if seen[idx]:
            continue
        tRes = [idx]
        for nbr in nbrLists[idx]:
            if not seen[nbr]:
                tRes.append(nbr)
                seen[nbr] = 1
        # update the number of neighbors:
        # remove all members of the new cluster from the list of
        # neighbors and reorder the tLists
        if reordering:
            # get the list of affected molecules, i.e. all molecules
            # which have at least one of the members of the new cluster
            # as a neighbor
            nbrNbr = [nbrLists[t] for t in tRes]
            nbrNbr = frozenset().union(*nbrNbr)
            # loop over all remaining molecules in tLists but only
            # consider unassigned and affected compounds
            for x, y in enumerate(tLists):
                y1 = y[1]
                if seen[y1] or (y1 not in nbrNbr):
                    continue
                # update the number of neighbors
                nbrLists[y1] = set(nbrLists[y1]).difference(tRes)
                tLists[x] = (len(nbrLists[y1]), y1)
            # now reorder the list
            tLists.sort(reverse=True)
        res.append(tuple(tRes))
    return tuple(res)

# Taylor butina using chemfp
def taylor_butina_cluster(arena, threshold):
    # Sort the results so that fingerprints with more hits come
    # first. This is more likely to be a cluster centroid. Break ties
    # arbitrarily by the fingerprint id; since fingerprints are
    # ordered by the number of bits this likely makes larger
    # structures appear first.:

    # Reorder so the centroid with the most hits comes first.  (That's why I do
    # a reverse search.)  Ignore the arbitrariness of breaking ties by
    # fingerprint index
    print "-------------------------------------------------"
    print "starting Butina clustering"
    # now cluster the data:
    start_time = time.time()

    similarity_table = search.threshold_tanimoto_search_symmetric(
        arena, threshold=threshold)

    centroid_table = sorted(((len(indices), i, indices)
                             for (i, indices) in enumerate(similarity_table.iter_indices())),
                            reverse=True)


    # Determine the true/false singletons and the clusters
    true_singletons = []
    false_singletons = []
    clusters = []

    seen = set()
    for (size, fp_idx, members) in centroid_table:
        if fp_idx in seen:
            # Can't use a centroid which is already assigned
            continue
        seen.add(fp_idx)

        # Figure out which ones haven't yet been assigned
        unassigned = set(members) - seen

        if not unassigned:
            mol_id = arena.ids[fp_idx]
            listi = [mol_id]
            false_singletons.append(listi)
            continue

        # this is a new cluster
        unassigned.add(fp_idx)

        #change indeces to ids
        unassigned_ids = []
        for mol in unassigned:
            unassigned_ids.append(arena.ids[mol])

        clusters.append(unassigned_ids)
        seen.update(unassigned)
    print "time taken: ", time.time() - start_time
    # Return the results:
    return ClusterResults(true_singletons, false_singletons, clusters)


#KMeans
def KMeansClustering(allFingerprints):
    print "-------------------------------------------------"
    print "starting KMeans clustering"
    arrayFingerprints = fingerprintsToNPArr(allFingerprints)

    start_time = time.time()
    kmeans = KMeans(n_clusters=10)
    kmeans.fit_predict(arrayFingerprints)
    print "time taken: ", time.time() - start_time
    labels = kmeans.labels_
    return labels

# Wards Clustering
# distance matrix, number of fingerprints
def WardsClustering(dists, nfps):
    print "-------------------------------------------------"
    print "starting Wards clustering"
    start_time = time.time()
    c_tree = Murtagh.ClusterData(dists, nfps, Murtagh.WARDS, isDistData=True)
    print "time taken: ", time.time() - start_time
    return c_tree


# get clusters at a particular hierarchy level
def GetHierarchicalLevel(cluster, levelNeeded, currentLevel=0):
    result_cluster = []
    # if this is the level required, or this level has no children, return the data points
    if ((currentLevel == levelNeeded) or (cluster.GetData() is not None)):
        points = cluster.GetPoints()
        for point in points:
            result_cluster.append(point.GetData())
        return result_cluster

    # recursive call for each child
    for child in cluster.GetChildren():
        result = GetHierarchicalLevel(child, levelNeeded, currentLevel + 1)
        # if a list of datapoints is returned, (ie the cluster) than append it
        if isinstance(result[0], int):
            result_cluster.append(result)
        else:
            # if list of clusters are returned, then extend list
            result_cluster.extend(result)
    return result_cluster
