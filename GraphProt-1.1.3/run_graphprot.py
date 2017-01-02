import os
import sys
import gzip
import pdb
from sklearn.metrics import roc_curve, auc, roc_auc_score
import subprocess
from shutil import copyfile

#copyfile(src, dst)

def load_label_seq(seq_file, outputfile, score_file, train = False): 
    label_list = []
    seq = ''
    if train:
    	fw = open(outputfile + 'posi.fa', 'w')
        fw1 = open(outputfile + 'nega.fa', 'w')
    else:
        fw = open(outputfile, 'w')
        
    fw2 = open(score_file, 'w')
    with gzip.open(seq_file, 'r') as fp:
        for line in fp:
            if line[0] == '>':
                name = line[1:-1]
                posi_label = name.split(';')[-1]
                posi_label1 = name.split(';')[0].split()[-1]
                label = posi_label.split(':')[-1]
                #pdb.set_trace()
                if int(label) < 0:
                    label = '0'
                if train:
                    if int(label):
                        fw.write('>'+ posi_label1 + '\n')
                    else:
                        fw1.write('>'+ posi_label1 + '\n')
                else: 
                    fw.write('>'+ posi_label1 + '\n')
                #label_list.append(int(label))
                fw2.write(label + '\n')
            else:
                seq = line[:-1].upper()
                if train:
                    if int(label):
                        fw.write(seq + '\n')
                    else:
                        fw1.write(seq + '\n')
                else:
                    fw.write(seq + '\n')
                
    fw.close()
    if train:
    	fw1.close()
    fw2.close()

def run_graphprot(protein, fw  = None):
    path = './datasets/clip/' + protein + '/30000/training_sample_0'
    load_label_seq(os.path.join(path, 'sequences.fa.gz'), 'train', 'scorefile', train = True)
    #pdb.set_trace()
    clistr = 'perl GraphProt.pl --action train -fasta trainposi.fa  -negfasta trainnega.fa'
    print clistr
    #fp = os.popen(clistr, 'r')
    #fp.close()
    score_str =subprocess.call(clistr, shell=True)
    #pdb.set_trace()
    path = './datasets/clip/' + protein + '/30000/test_sample_0'
    load_label_seq(os.path.join(path, 'sequences.fa.gz'), 'test.fa', 'scoretest')
    
    clistr = 'perl GraphProt.pl --action predict -fasta test.fa -model GraphProt.model'
    #fp = os.popen(clistr, 'r')
    #fp.close()
    score_str =subprocess.call(clistr, shell=True)
    #pdb.set_trace()
    true_y = []
    with open('scoretest' , 'r') as fp:
        for line in fp:
            values = line.rstrip().split()
            true_y.append(int(values[0]))
    predict_prob = []        
    with open('GraphProt.predictions', 'r') as fp:
        for line in fp:
            values = line.rstrip().split()
            predict_prob.append(float(values[-1]))
                    
    auc = roc_auc_score(true_y, predict_prob)
    
    print "Test AUC: ", auc    
    fw.write(str(auc) + '\n')
    mylabel = "\t".join(map(str, true_y))
    myprob = "\t".join(map(str, predict_prob))
    fw.write(mylabel + '\n')
    fw.write(myprob + '\n')
    copyfile('GraphProt.predictions', protein + '_GraphProt.predictions')

def run_predict():
    data_dir = './datasets/clip'
    fw = open('result_classification_graphport', 'w')
    for protein in os.listdir(data_dir):
        print protein
        fw.write(protein + '\t')
        run_graphprot(protein, fw= fw)
    fw.close()
    
run_predict()
