from Clustering import *
from Evaluation import *
from Reader import *

QPI = True
F1 = True

if __name__ == "__main__":
    actives_file = '../mols/compounds5.smi'
    results_file = '../mols/resultsSerial/resultsclsmi.smi'

    actives = read_smiles_fingerprints(actives_file)
    clusters = read_smiles_fps_clusters(results_file)
    #for cluster in clusters:
            #print len(cluster)

    if F1:
        print ('---------------------------------------F Score---------------------------------------')

        max_score, all_scores = calculate_f1(clusters, actives)
        for score in all_scores:
            if(score[3] != 0):
                print ('Cluster %s: Precision %s Recall %s FScore %s Having %s actives %d no activity mols' % (score[0],score[1],score[2],score[3],score[4],score[5]))

        print ('Max Result : Precision %s Recall %s FScore %s' % (max_score[0], max_score[1], max_score[2]))

    if QPI:
        print ('---------------------------------------QPI---------------------------------------')