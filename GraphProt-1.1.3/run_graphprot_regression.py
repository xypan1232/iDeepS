import os
import sys
import gzip
import pdb
from sklearn.metrics import roc_curve, auc, roc_auc_score
import subprocess

def load_label_seq(seq_file, outputfile, score_file): 
    label_list = []
    seq = ''
    fw = open(outputfile, 'w')
    fw2 = open(score_file, 'w')
    with gzip.open(seq_file, 'r') as fp:
        for line in fp:
            if line[0] == '>':
                name = line[1:-1]
                #fw.write(line)
                posi_label = name.split(';')[-1]
                fw.write('>' + name.strip().split(';')[0] + '\n')
                #pdb.set_trace()
                label = posi_label.split(':')[-1]
                if int(label) < 0:
                    label = '0'
                #label_list.append(int(label))
                fw2.write(label + '\n')
            else:
                seq = line[:-1].upper()
                fw.write(seq + '\n')
    fw.close()
    fw2.close()

def run_graphprot(protein, fw  = None):
    path = './datasets/clip/' + protein + '/30000/training_sample_0'
    load_label_seq(os.path.join(path, 'sequences.fa.gz'), 'train.fa', 'scorefile')
    #pdb.set_trace()
    clistr = 'perl GraphProt.pl --action train -mode regression -fasta train.fa -affinities scorefile'
    print clistr
    score_str =subprocess.call(clistr, shell=True)
    #fp = os.popen(clistr, 'r')
    #fp.close()
    #pdb.set_trace()
    path = './datasets/clip/' + protein + '/30000/test_sample_0'
    load_label_seq(os.path.join(path, 'sequences.fa.gz'), 'test.fa', 'scoretest')
    
    clistr = 'perl GraphProt.pl --action predict -mode regression -fasta test.fa -model GraphProt.model'
    
    score_str =subprocess.call(clistr, shell=True)
    true_y = []
    with open('scoretest' , 'r') as fp:
        for line in fp:
            values = line.rstrip().split()
            true_y.append(int(values[0]))
    predict_prob = []        
    with open('GraphProt.predictions', 'r') as fp:
        for line in fp:
            values = line.rstrip().split()
            predict_prob.append(float(values[1]))
                    
    auc = roc_auc_score(true_y, predict_prob)
    
    print "Test AUC: ", auc    
    fw.write(str(auc) + '\n')

def run_predict():
    data_dir = './datasets/clip'
    fw = open('result_file_graphport', 'w')
    for protein in os.listdir(data_dir):
        print protein
        fw.write(protein + '\t')
        run_graphprot(protein, fw= fw)
    fw.close()
    
run_predict()
