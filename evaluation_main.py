from Clustering import *
from Evaluation import *
from Reader import *
from os import listdir
from os.path import isfile,join

QPI = True
F1 = True
EF = True

def calculate_f1_results(clusters, actives, show_all= True):
    print ('---------------------------------------F Score---------------------------------------')

    max_score, all_scores = calculate_f1(clusters, actives)

    if show_all:
        for score in all_scores:
            if (score[3] != 0):
                print ('Cluster %s: Precision %s Recall %s FScore %s Having %s actives %d no activity mols' % (
                score[0], score[1], score[2], score[3], score[4], score[5]))

    print ('Max Result : Precision %s Recall %s FScore %s' % (max_score[0], max_score[1], max_score[2]))

def calculate_qpi_results(clusters, actives):
    print ('---------------------------------------QPI---------------------------------------')
    dataset_size = 0
    for cluster in clusters:
        dataset_size += len(cluster)

    qpi = calculate_qpi(clusters, actives, dataset_size, 'RDKIT')
    print "QPI:", qpi

def calculate_ef_results(clusters, actives, show_all= True):
    print ('---------------------------------------Enrichment Factor---------------------------------------')

    max_score, all_scores = calculate_enrichment(clusters, actives)
    if show_all:
        for score in all_scores:
            if (score[1] != 0):
                print ('Cluster %s: Enrichment %s Having %s actives %d no activity mols' % (
                score[0], score[1], score[2], score[3]))

    print ('Max Result : Enrichment %s' % (max_score))

def evaluate_all_files(actives_file, directory_path, show_all=True):
    onlyfiles = [join(directory_path, f) for f in listdir(directory_path) if isfile(join(directory_path, f))]

    actives = read_smiles_fingerprints(actives_file)
    for file in onlyfiles:
        clusters = read_smiles_fps_clusters(file)
        print (file)
        if F1:
            calculate_f1_results(clusters, actives,show_all)

        if QPI:
            calculate_qpi_results(clusters, actives)

        if EF:
            calculate_ef_results(clusters, actives,show_all)

def evaluate_single_file(actives_file, results_file):
    actives = read_smiles_fingerprints(actives_file)
    clusters = read_smiles_fps_clusters(results_file)

    if F1:
        calculate_f1_results(clusters, actives)

    if QPI:
        calculate_qpi_results(clusters, actives)

    if EF:
        calculate_ef_results(clusters, actives)

if __name__ == "__main__":
    actives_file = '../mols/SmilesActives/ABL1.smi'
    results_file = '../mols/resultsSerial/butina/butinaResult0.3.smi'
    directory_path = '../mols/resultsSerial/wards'
    #results_file = '../mols/resultsSpark/result/part-00000'

    evaluate_all_files(actives_file, directory_path, False)

    #actives = read_smiles_fingerprints(actives_file)
    #clusters = read_smiles_fps_clusters(results_file)

    #if F1:
    #    calculate_f1(clusters, actives)

    #if QPI:
    #    calculate_qpi(clusters, actives)

    #if EF:
    #   calculate_ef(clusters, actives)
