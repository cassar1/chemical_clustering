# adapted from http://www.dalkescientific.com/writings/taylor_butina.py
from rdkit import DataStructs
import time
from Reader import *
import sys
import copy
import numpy
from chemfp import search
import random

# The results of the Taylor-Butina clustering
class ClusterResults(object):
    def __init__(self, true_singletons, false_singletons, clusters):
        self.true_singletons = true_singletons
        self.false_singletons = false_singletons
        self.clusters = clusters


def combine_clusters(cluster_results):
    clusters = cluster_results.clusters

    true_singletons = cluster_results.true_singletons
    false_singletons = cluster_results.false_singletons
    clusters.extend(true_singletons)
    clusters.extend(false_singletons)

    return clusters


def distance_matrix(arena):
    n = len(arena)

    # Start off a similarity matrix with 1.0s along the diagonal
    similarities = numpy.identity(n, "d")

    ## Compute the full similarity matrix.
    # The implementation computes the upper-triangle then copies
    # the upper-triangle into lower-triangle. It does not include
    # terms for the diagonal.
    results = search.threshold_tanimoto_search_symmetric(arena, threshold=0.0)

    # Copy the results into the NumPy array.
    for row_index, row in enumerate(results.iter_indices_and_scores()):
        for target_index, target_score in row:
            similarities[row_index, target_index] = target_score

    # Return the distance matrix using the similarity matrix
    return 1.0 - similarities


# output cluster results object
def report_cluster_results(cluster_results, arena):
    true_singletons = cluster_results.true_singletons
    false_singletons = cluster_results.false_singletons
    clusters = cluster_results.clusters

    # Sort the singletons by id.
    print(len(true_singletons), "true singletons")
    print("=>", " ".join(sorted(arena.ids[idx] for idx in true_singletons)))

    print(len(false_singletons), "false singletons")

    # print("=>", " ".join(sorted(arena.ids[idx] for idx in false_singletons)))

    # Sort so the cluster with the most compounds comes first,
    # then by alphabetically smallest id
    def cluster_sort_key(cluster):
        centroid_id, cluster_members = cluster
        return -len(cluster_members), arena.ids[centroid_id]

    clusters.sort(key=cluster_sort_key)

    total_compounds = 0
    print(len(clusters), "clusters")
    for centroid_idx, members in clusters:
        print(arena.ids[centroid_idx], "has", len(members), "other members")
        total_compounds = total_compounds + len(members) + 1
        # print("=>", " ".join(arena.ids[idx] for idx in members))
    print "Total Compounds = ", total_compounds + len(false_singletons)


# given a list of active and decoys, returns a percentage of the list
def getPercentageMolecules(decimal, activeDescriptors, decoyDescriptors, shuffle=False):
    activeAmount = len(activeDescriptors)
    decoyAmount = len(decoyDescriptors)

    activeAmount = int(activeAmount * decimal)
    decoyAmount = int(decoyAmount* decimal)

    if shuffle:
        random.shuffle(activeDescriptors)
        random.shuffle(decoyDescriptors)
    activeDescriptorsSubset = activeDescriptors[0:activeAmount]
    decoyDescriptorsSubset = decoyDescriptors[0:decoyAmount]

    return activeDescriptorsSubset, decoyDescriptorsSubset


# given an active arena and arena of all molecules, returns a percentage of the arena
def get_arena_percentage(decimal, arena_active, arena_all, shuffle=False):
    active_amount = len(arena_active)
    total_decoy_amount = len(arena_all) - len(arena_active)

    active_amount = int(active_amount * decimal)
    decoy_amount = int(total_decoy_amount * decimal)

    # iterate active arena, and for each molecule, get its relative index in the arena containing all molecules
    active_indeces = []
    for item in arena_active:
        # index = arena_all.get_index_by_id(item)
        active_mols = search.contains_fp(item[1], arena_all)
        for act in active_mols:
            if act[0] not in active_indeces:
                active_indeces.append(act[0])

    decoy_counter = 0
    added_decoys_amount = 0
    decoy_indeces = []
    while added_decoys_amount < total_decoy_amount:
        if decoy_counter not in active_indeces:
            decoy_indeces.append(decoy_counter)
            added_decoys_amount = added_decoys_amount + 1
        decoy_counter = decoy_counter + 1

    if shuffle:
        random.shuffle(active_indeces)
        random.shuffle(decoy_indeces)

    active_indeces_subset = active_indeces[0:active_amount]

    molecule_indeces = active_indeces_subset[:]
    molecule_indeces.extend(decoy_indeces[0:decoy_amount])

    arena_subset = arena_all.copy(molecule_indeces)
    # print "Active subset list size ", len(active_indeces_subset)
    arena_active_subset = arena_all.copy(active_indeces_subset)

    return arena_active_subset, arena_subset

#return similarities as one d list
def CalcSimilarities(fps):
    #first generate the distance matrix:
    print "starting distance calculations"
    start_time = time.time()
    dists = []
    dists_2d = []
    nfps = len(fps)
    print "total fingerprints ", nfps
    print "Type of : ", type(0.56)
    print "Size of double: ", sys.getsizeof(0.56)
    print "Size of list double: ", sys.getsizeof([0.56])
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        #print "length result ", len(sims) , " memory size ", sys.getsizeof(sims) , " bytes "
        dists.extend([1-x for x in sims])
        dists_2d.append(sims)
        #print "Similarity store size: ", sys.getsizeof(dists)
    print "time taken: ", time.time() - start_time
    return dists, dists_2d


# calculate neighbour list using RDKIT
#adapted from RDKIT BUTINA Algorithm
def get_neighbours_list(fps, threshold):
    start_time = time.time()
    n_points = len(fps)
    neighbour_list = [None] * n_points
    for i in range(n_points):
        neighbour_list[i] = []

    for i in range(n_points):
        #print "chemical", i
        for j in range(i):
            if i != j:
                distance = DataStructs.FingerprintSimilarity(fps[i], fps[j], metric=DataStructs.TanimotoSimilarity)

                if (distance >= threshold):
                    # print "comparing to chemical", j , "distance ", distance
                    neighbour_list[i].append(j)
                    neighbour_list[j].append(i)
    # sort list of neighbours by num neighbours
    tLists = [(len(y), x) for x, y in enumerate(neighbour_list)]
    tLists.sort(reverse=True)
    print "time taken to calculate ", n_points ," : ", time.time() - start_time
    return tLists, neighbour_list

# Just a test method to output the bits of the molecules and a count of the duplicates
def output_bit_strings(input_file):
    bitstring = ""
    oldbitstring = ""
    numberDuplicated = 0
    activeDescriptors = read_fingerprints(input_file)

    for active in activeDescriptors:
        bitstring = ""
        for i in range(0, 1024):
            bitstring += str(active[i])
        print bitstring

        if bitstring == oldbitstring:
            print "duplicate found"
            numberDuplicated = numberDuplicated + 1
        oldbitstring = bitstring
    print "Total actives ", len(activeDescriptors)
    print "non duplicates = ", len(activeDescriptors) - numberDuplicated
    print numberDuplicated, " Duplicated Found"

def distance_matrix(arena):
    print "Here 1"
    n = len(arena)
    # Start off a similarity matrix with 1.0s along the diagonal
    similarities = numpy.identity(n, "d")

    # Compute the full similarity matrix.
    # The implementation computes the upper-triangle then copies
    # the upper-triangle into lower-triangle. It does not include
    # terms for the diagonal.
    results = search.threshold_tanimoto_search_symmetric(arena, threshold=0.0,include_lower_triangle=False)

    # Copy the results into the NumPy array.
    for row_index, row in enumerate(results.iter_indices_and_scores()):
        for target_index, target_score in row:
            similarities[row_index, target_index] = target_score

    # Return the distance matrix using the similarity matrix
    return 1.0 - similarities