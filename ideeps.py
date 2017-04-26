import sys
import os
import numpy
import pdb
from keras.models import Sequential, model_from_config
from keras.layers.core import Dense, Dropout, Activation, Flatten, Merge
#from keras.layers import Input, merge, LSTM
from keras.layers.normalization import BatchNormalization
from keras.layers.advanced_activations import PReLU
from keras.utils import np_utils, generic_utils
from keras.optimizers import SGD, RMSprop, Adadelta, Adagrad, Adam
from keras.layers import normalization
from keras.layers.convolutional import Convolution2D, MaxPooling2D
from keras.layers import LSTM, Bidirectional 
from keras.layers.embeddings import Embedding
from keras.layers.convolutional import Convolution2D, MaxPooling2D,Convolution1D, MaxPooling1D
from keras import regularizers
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras.constraints import maxnorm
#from seya.layers.recurrent import Bidirectional
from sklearn import svm, grid_search
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.cross_validation import train_test_split
from sklearn.calibration import CalibratedClassifierCV
from sklearn.cross_validation import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import random
import gzip
from sklearn import svm
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import PCA
from sklearn import metrics
from sklearn.metrics import roc_auc_score
from sklearn.cross_validation import train_test_split
from sklearn.grid_search import GridSearchCV
from scipy import sparse
import pdb
from math import  sqrt
from sklearn.metrics import roc_curve, auc
import theano
import subprocess as sp
import scipy.stats as stats
#from seq_motifs import *
#import structure_motifs
from keras import backend as K
from rnashape_structure import run_rnashape

def calculate_performace(test_num, pred_y,  labels):
    tp =0
    fp = 0
    tn = 0
    fn = 0
    for index in range(test_num):
        if labels[index] ==1:
            if labels[index] == pred_y[index]:
                tp = tp +1
            else:
                fn = fn + 1
        else:
            if labels[index] == pred_y[index]:
                tn = tn +1
            else:
                fp = fp + 1               
            
    acc = float(tp + tn)/test_num
    precision = float(tp)/(tp+ fp)
    sensitivity = float(tp)/ (tp+fn)
    specificity = float(tn)/(tn + fp)
    MCC = float(tp*tn-fp*fn)/(np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return acc, precision, sensitivity, specificity, MCC 

def merge_seperate_network(X_train1, X_train2, Y_train):
    left_hid = 128
    right_hid = 64
    left = get_rnn_fea(X_train1, sec_num_hidden = left_hid)
    right = get_rnn_fea(X_train2, sec_num_hidden = right_hid)
    
    model = Sequential()
    model.add(Merge([left, right], mode='concat'))
    total_hid = left_hid + right_hid
    
    model.add(Dense(total_hid, 2))
    model.add(Dropout(0.3))
    model.add(Activation('softmax'))
    
    sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
    model.compile(loss='categorical_crossentropy', optimizer=sgd) #'rmsprop')
    
    model.fit([X_train1, X_train2], Y_train, batch_size=100, nb_epoch=100, verbose=0)
    
    return model

def read_seq(seq_file):
    seq_list = []
    seq = ''
    with gzip.open(seq_file, 'r') as fp:
        for line in fp:
            if line[0] == '>':
                name = line[1:-1]
                if len(seq):
                    seq_array = get_RNA_seq_concolutional_array(seq)
                    seq_list.append(seq_array)                    
                seq = ''
            else:
                seq = seq + line[:-1]
        if len(seq):
            seq_array = get_RNA_seq_concolutional_array(seq)
            seq_list.append(seq_array) 
    
    return np.array(seq_list)

def load_label_seq(seq_file):
    label_list = []
    seq = ''
    with gzip.open(seq_file, 'r') as fp:
        for line in fp:
            if line[0] == '>':
                name = line[1:-1]
                posi_label = name.split(';')[-1]
                label = posi_label.split(':')[-1]
                label_list.append(int(label))
    return np.array(label_list)

def read_rnashape(structure_file):
    struct_dict = {}
    with gzip.open(structure_file, 'r') as fp:
        for line in fp:
            if line[0] == '>':
                name = line[:-1]
            else:
                strucure = line[:-1]
                struct_dict[name] = strucure
                
    return struct_dict

def read_structure(seq_file, path):
    seq_list = []
    structure_list = []
    struct_exist = False
    if not os.path.exists(path + '/structure.gz'):
        fw = gzip.open(path + '/structure.gz', 'w')
    else:
        fw = None
        struct_exist = True
        struct_dict = read_rnashape(path + '/structure.gz')
        #pdb.set_trace()
    seq = ''
    with gzip.open(seq_file, 'r') as fp:
        for line in fp:
            if line[0] == '>':
                name = line
                if len(seq):
                    if struct_exist:
                        structure = struct_dict[old_name[:-1]]
                        seq_array, struct = get_RNA_structure_concolutional_array(seq, fw, structure = structure)
                    else:
                        fw.write(old_name)
                        seq_array, struct = get_RNA_structure_concolutional_array(seq, fw)
                    seq_list.append(seq_array)
                    structure_list.append(struct)
                old_name = name              
                seq = ''
            else:
                seq = seq + line[:-1]
        if len(seq): 
            if struct_exist:
                structure = struct_dict[old_name[:-1]]
                seq_array, struct = get_RNA_structure_concolutional_array(seq, fw, structure = structure)
            else:
                fw.write(old_name)
                seq_array, struct = get_RNA_structure_concolutional_array(seq, fw)
            #seq_array, struct = get_RNA_structure_concolutional_array(seq, fw)
            seq_list.append(seq_array)
            structure_list.append(struct)  
    if fw:
        fw.close()
    return np.array(seq_list), structure_list


def read_oli_feature(seq_file):
    trids4 = get_4_trids()
    seq_list = []
    seq = ''
    with gzip.open(seq_file, 'r') as fp:
        for line in fp:
            if line[0] == '>':
                name = line[1:-1]
                if len(seq):
                    seq_array = get_4_nucleotide_composition(trids4, seq)
                    seq_list.append(seq_array)                    
                seq = ''
            else:
                seq = seq + line[:-1]
        if len(seq):
            seq_array = get_4_nucleotide_composition(trids4, seq)
            seq_list.append(seq_array) 
    
    return np.array(seq_list)    

def get_4_trids():
    nucle_com = []
    chars = ['A', 'C', 'G', 'U']
    base=len(chars)
    end=len(chars)**4
    for i in range(0,end):
        n=i
        ch0=chars[n%base]
        n=n/base
        ch1=chars[n%base]
        n=n/base
        ch2=chars[n%base]
        n=n/base
        ch3=chars[n%base]
        nucle_com.append(ch0 + ch1 + ch2 + ch3)
    return  nucle_com

def get_4_nucleotide_composition(tris, seq, pythoncount = True):
    #pdb.set_trace()
    seq_len = len(seq)
    seq = seq.upper().replace('T', 'U')
    tri_feature = []
    
    if pythoncount:
        for val in tris:
            num = seq.count(val)
            tri_feature.append(float(num)/seq_len)
    else:
        k = len(tris[0])
        tmp_fea = [0] * len(tris)
        for x in range(len(seq) + 1- k):
            kmer = seq[x:x+k]
            if kmer in tris:
                ind = tris.index(kmer)
                tmp_fea[ind] = tmp_fea[ind] + 1
        tri_feature = [float(val)/seq_len for val in tmp_fea]
        #pdb.set_trace()        
    return tri_feature

def load_data(path, seq = True, oli = True):
    """
        Load data matrices from the specified folder.
    """

    data = dict()
    if seq: 
        tmp = []
        tmp.append(read_seq(os.path.join(path, 'sequences.fa.gz')))
        seq_onehot, structure = read_structure(os.path.join(path, 'sequences.fa.gz'), path)
        tmp.append(seq_onehot)
        data["seq"] = tmp
        data["structure"] = structure
    
    if oli: data["oli"] = read_oli_feature(os.path.join(path, 'sequences.fa.gz'))
    
    data["Y"] = load_label_seq(os.path.join(path, 'sequences.fa.gz'))
    #np.loadtxt(gzip.open(os.path.join(path,
                #                            "matrix_Response.tab.gz")),
                #                            skiprows=1)
    #data["Y"] = data["Y"].reshape((len(data["Y"]), 1))

    return data   

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    complseq = [complement[base] for base in seq]
    return complseq

def reverse_complement(seq):
    seq = list(seq)
    seq.reverse()
    return ''.join(complement(seq))

def preprocess_data(X, scaler=None, stand = False):
    if not scaler:
        if stand:
            scaler = StandardScaler()
        else:
            scaler = MinMaxScaler()
        scaler.fit(X)
    X = scaler.transform(X)
    return X, scaler    

def preprocess_labels(labels, encoder=None, categorical=True):
    if not encoder:
        encoder = LabelEncoder()
        encoder.fit(labels)
    y = encoder.transform(labels).astype(np.int32)
    if categorical:
        y = np_utils.to_categorical(y)
    return y, encoder

def get_RNA_seq_concolutional_array(seq, motif_len = 4):
    seq = seq.replace('U', 'T')
    alpha = 'ACGT'
    #for seq in seqs:
    #for key, seq in seqs.iteritems():
    row = (len(seq) + 2*motif_len - 2)
    new_array = np.zeros((row, 4))
    for i in range(motif_len-1):
        new_array[i] = np.array([0.25]*4)
    
    for i in range(row-3, row):
        new_array[i] = np.array([0.25]*4)
        
    #pdb.set_trace()
    for i, val in enumerate(seq):
        i = i + motif_len-1
        if val not in 'ACGT':
            new_array[i] = np.array([0.25]*4)
            continue
        #if val == 'N' or i < motif_len or i > len(seq) - motif_len:
        #    new_array[i] = np.array([0.25]*4)
        #else:
        try:
            index = alpha.index(val)
            new_array[i][index] = 1
        except:
            pdb.set_trace()
        #data[key] = new_array
    return new_array

def get_RNA_structure_concolutional_array(seq, fw, structure = None, motif_len = 6):
    if fw is None:
        struc_en = structure
    else:
        #print 'running rnashapes'
        seq = seq.replace('U', 'T')
        struc_en = run_rnashape(seq)
        fw.write(struc_en + '\n')
        
    alpha = 'FTIHMS'
    row = (len(struc_en) + 2*motif_len - 2)
    new_array = np.zeros((row, 6))
    for i in range(motif_len-1):
        new_array[i] = np.array([0.16]*6)
    
    for i in range(row-5, row):
        new_array[i] = np.array([0.16]*6)

    for i, val in enumerate(struc_en):
        i = i + motif_len-1
        if val not in alpha:
            new_array[i] = np.array([0.16]*6)
            continue
        try:
            index = alpha.index(val)
            new_array[i][index] = 1
        except:
            pdb.set_trace()
        
    return new_array, struc_en

def get_2d_cnn_network():
    nb_conv = 4
    nb_pool = 2
    model = Sequential()
    model.add(Convolution2D(64, nb_conv, nb_conv,
                            border_mode='valid',
                            input_shape=(1, 107, 4)))
    model.add(Activation('relu'))
    model.add(Convolution2D(64, nb_conv, nb_conv))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(nb_pool, nb_pool)))
    model.add(Dropout(0.25))
    model.add(Flatten())
    
    return model

def fork (model, n=2):
    forks = []
    for i in range(n):
        f = Sequential()
        f.add (model)
        forks.append(f)
    return forks

'''
class MyReshape(Layer):
    def get_output(self, train):
        X = self.get_input(train)
        nshape = (1,) + X.shape 
        return theano.tensor.reshape(X, nshape)
'''

def set_cnn_model(input_dim, input_length):
    nbfilter = 16
    model = Sequential()
    model.add(Convolution1D(input_dim=input_dim,input_length=input_length,
                            nb_filter=nbfilter,
                            filter_length=10,
                            border_mode="valid",
                            #activation="relu",
                            subsample_length=1))
    model.add(Activation('relu'))
    model.add(MaxPooling1D(pool_length=3))

    model.add(Dropout(0.5))

    return model

def get_cnn_network():
    '''
     get_feature = theano.function([origin_model.layers[0].input],origin_model.layers[11].get_output(train=False),allow_input_downcast=False)
    feature = get_feature(data)
    '''
    nbfilter = 16
    print 'configure cnn network'

    seq_model = set_cnn_model(4, 107)
    struct_model = set_cnn_model(6, 111)
    #pdb.set_trace()
    model = Sequential()
    model.add(Merge([seq_model, struct_model], mode='concat', concat_axis=1))

    model.add(Bidirectional(LSTM(2*nbfilter)))
    
    model.add(Dropout(0.10))
    
    model.add(Dense(nbfilter*2, activation='relu'))
    print model.output_shape
    
    return model

def get_cnn_network_old():
    '''
     get_feature = theano.function([origin_model.layers[0].input],origin_model.layers[11].get_output(train=False),allow_input_downcast=False)
    feature = get_feature(data)
    '''
    print 'configure cnn network'
    nbfilter = 16


    model = Sequential()
    model.add(Convolution1D(input_dim=4,input_length=107,
                            nb_filter=nbfilter,
                            filter_length=10,
                            border_mode="valid",
                            #activation="relu",
                            subsample_length=1))
    model.add(Activation('relu'))
    model.add(MaxPooling1D(pool_length=3))
    
    model.add(Dropout(0.5))
    
    model.add(Flatten())
    
    model.add(Dense(nbfilter, activation='relu'))

    model.add(Dropout(0.25))

    
    return model

def get_struct_network():
    '''
     get_feature = theano.function([origin_model.layers[0].input],origin_model.layers[11].get_output(train=False),allow_input_downcast=False)
    feature = get_feature(data)
    '''
    print 'configure cnn network'
    nbfilter = 16

    model = Sequential()
    model.add(Convolution1D(input_dim=6,input_length=107,
                            nb_filter=nbfilter,
                            filter_length=10,
                            border_mode="valid",
                            #activation="relu",
                            subsample_length=1))
    model.add(Activation('relu'))
    model.add(MaxPooling1D(pool_length=3))
    
    model.add(Dropout(0.5))

    
    model.add(Flatten())
    
    model.add(Dense(nbfilter, activation='relu'))

    model.add(Dropout(0.25))
    
    return model


def get_rnn_fea(train, sec_num_hidden = 128, num_hidden = 128):
    print 'configure network for', train.shape
    model = Sequential()

    model.add(Dense(num_hidden, input_shape=(train.shape[1],), activation='relu'))
    model.add(PReLU())
    model.add(BatchNormalization(mode=2))
    model.add(Dropout(0.5))
    model.add(Dense(num_hidden, input_dim=num_hidden, activation='relu'))
    #model.add(Dense(num_hidden, input_shape=(num_hidden,), activation='relu'))
    model.add(PReLU())
    model.add(BatchNormalization(mode=2))
    #model.add(Activation('relu'))
    model.add(Dropout(0.5))

    return model

def get_structure_motif_fig(filter_weights, filter_outs, out_dir, protein, seq_targets, sample_i = 0, structure = None):
    print 'plot motif fig', out_dir
    #seqs, seq_targets = get_seq_targets(protein)
    seqs = structure
    if sample_i:
        print 'sampling'
        seqs = []
        for ind, val in enumerate(seqs):
            if ind in sample_i:
                seqs.append(val)
            
        
        seq_targets = seq_targets[sample_i]
        filter_outs = filter_outs[sample_i]
    
    num_filters = filter_weights.shape[0]
    filter_size = 7 #filter_weights.shape[2]

    
    #################################################################
    # individual filter plots
    #################################################################
    # also save information contents
    filters_ic = []
    meme_out = structure_motifs.meme_intro('%s/filters_meme.txt'%out_dir, seqs)

    for f in range(num_filters):
        print 'Filter %d' % f

        # plot filter parameters as a heatmap
        structure_motifs.plot_filter_heat(filter_weights[f,:,:], '%s/filter%d_heat.pdf' % (out_dir,f))

        # write possum motif file
        structure_motifs.filter_possum(filter_weights[f,:,:], 'filter%d'%f, '%s/filter%d_possum.txt'%(out_dir,f), False)
        
        structure_motifs.plot_filter_logo(filter_outs[:,:, f], filter_size, seqs, '%s/filter%d_logo'%(out_dir,f), maxpct_t=0.5)
        
        filter_pwm, nsites = structure_motifs.make_filter_pwm('%s/filter%d_logo.fa'%(out_dir,f))
        if nsites < 10:
            # no information
            filters_ic.append(0)
        else:
            # compute and save information content
            filters_ic.append(info_content(filter_pwm))

            # add to the meme motif file
            structure_motifs.meme_add(meme_out, f, filter_pwm, nsites, False)

    meme_out.close()
    
            
def get_motif_fig(filter_weights, filter_outs, out_dir, protein, sample_i = 0):
    print 'plot motif fig', out_dir
    seqs, seq_targets = get_seq_targets(protein)
    if sample_i:
        print 'sampling'
        seqs = []
        for ind, val in enumerate(seqs):
            if ind in sample_i:
                seqs.append(val)
            
        
        seq_targets = seq_targets[sample_i]
        filter_outs = filter_outs[sample_i]
    
    num_filters = filter_weights.shape[0]
    filter_size = 7#filter_weights.shape[2]

    #pdb.set_trace()
    #################################################################
    # individual filter plots
    #################################################################
    # also save information contents
    filters_ic = []
    meme_out = meme_intro('%s/filters_meme.txt'%out_dir, seqs)

    for f in range(num_filters):
        print 'Filter %d' % f

        # plot filter parameters as a heatmap
        plot_filter_heat(filter_weights[f,:,:], '%s/filter%d_heat.pdf' % (out_dir,f))

        # write possum motif file
        filter_possum(filter_weights[f,:,:], 'filter%d'%f, '%s/filter%d_possum.txt'%(out_dir,f), False)

        # plot weblogo of high scoring outputs
        plot_filter_logo(filter_outs[:,:, f], filter_size, seqs, '%s/filter%d_logo'%(out_dir,f), maxpct_t=0.5)

        # make a PWM for the filter
        filter_pwm, nsites = make_filter_pwm('%s/filter%d_logo.fa'%(out_dir,f))

        if nsites < 10:
            # no information
            filters_ic.append(0)
        else:
            # compute and save information content
            filters_ic.append(info_content(filter_pwm))

            # add to the meme motif file
            meme_add(meme_out, f, filter_pwm, nsites, False)

    meme_out.close()


    #################################################################
    # annotate filters
    #################################################################
    # run tomtom #-evalue 0.01 
    subprocess.call('tomtom -dist pearson -thresh 0.05 -eps -oc %s/tomtom %s/filters_meme.txt %s' % (out_dir, out_dir, 'Ray2013_rbp_RNA.meme'), shell=True)

    # read in annotations
    filter_names = name_filters(num_filters, '%s/tomtom/tomtom.txt'%out_dir, 'Ray2013_rbp_RNA.meme')


    #################################################################
    # print a table of information
    #################################################################
    table_out = open('%s/table.txt'%out_dir, 'w')

    # print header for later panda reading
    header_cols = ('', 'consensus', 'annotation', 'ic', 'mean', 'std')
    print >> table_out, '%3s  %19s  %10s  %5s  %6s  %6s' % header_cols

    for f in range(num_filters):
        # collapse to a consensus motif
        consensus = filter_motif(filter_weights[f,:,:])

        # grab annotation
        annotation = '.'
        name_pieces = filter_names[f].split('_')
        if len(name_pieces) > 1:
            annotation = name_pieces[1]

        # plot density of filter output scores
        fmean, fstd = plot_score_density(np.ravel(filter_outs[:,:, f]), '%s/filter%d_dens.pdf' % (out_dir,f))

        row_cols = (f, consensus, annotation, filters_ic[f], fmean, fstd)
        print >> table_out, '%-3d  %19s  %10s  %5.2f  %6.4f  %6.4f' % row_cols

    table_out.close()


    #################################################################
    # global filter plots
    #################################################################
    if True:
        new_outs = []
        for val in filter_outs:
            new_outs.append(val.T)
        filter_outs = np.array(new_outs)
        print filter_outs.shape
        # plot filter-sequence heatmap
        plot_filter_seq_heat(filter_outs, '%s/filter_seqs.pdf'%out_dir)

def calculateInformationContentOfMostLikelySequence(self, sequenceNumber):
        mostLikelyPath = self.getMostLikelyPath()[1:-1]

        informationContent = []
        for stateIndex in mostLikelyPath:
            distribution = self.getEmission(stateIndex)

            shannonEntropy = 0
            for baseIndex in xrange(4):
                shannonEntropy -= distribution[baseIndex] * math.log(distribution[baseIndex], 2)
            smallSampleCorrection = 3.0 / (math.log1p(2) * 2 * sequenceNumber)
            informationContent.append(2 - (shannonEntropy + smallSampleCorrection))

        return informationContent

def calculateInformationContentStructure(self, sequenceNumber):
    probabilityToGetIntoState = self._getStateProbabilities()

    informationContent = []
    for position in xrange((self.N - 2 ) / 5):
        probs = []
        for structureIndex in xrange(5):
            probs.append(probabilityToGetIntoState[1 + (5 * position) + structureIndex])

        shannonEntropy = 0
        for structureIndex in xrange(5):
            scaledProb = probs[structureIndex] / float(sum(probs))
            shannonEntropy -= scaledProb * math.log(scaledProb, 2)
        smallSampleCorrection = 4.0 / (math.log1p(2) * 2 * sequenceNumber)
        informationContent.append(math.log(5, 2) - (shannonEntropy + smallSampleCorrection))

    return informationContent

def calculateInformationContentCombined(self, sequenceNumber):
    probabilityToGetIntoState = self._getStateProbabilities()

    informationContent = []
    for position in xrange((self.N - 2 ) / 5):
        probs = []
        for stateIndex in xrange(1 + (5 * position), 6 + (5 * position)):
            distribution = self.getEmission(stateIndex)
            for baseIndex in xrange(4):
                probs.append(probabilityToGetIntoState[stateIndex] * distribution[baseIndex])

        shannonEntropy = 0
        for combinedIndex in xrange(20):
            scaledProb = probs[combinedIndex] / float(sum(probs))
            shannonEntropy -= scaledProb * math.log(scaledProb, 2)
        smallSampleCorrection = 19.0 / (math.log1p(2) * 2 * sequenceNumber)
        informationContent.append(math.log(20, 2) - (shannonEntropy + smallSampleCorrection))

    return informationContent
    
def get_seq_targets(protein):
    path = "./datasets/clip/%s/30000/test_sample_0" % protein
    data = load_data(path)
    seq_targets = np.array(data['Y'])
    
    seqs = []
    seq = ''
    fp = gzip.open(path +'/sequences.fa.gz')
    for line in fp:
        if line[0] == '>':
            name = line[1:-1]
            if len(seq):
                seqs.append(seq)                    
            seq = ''
        else:
            seq = seq + line[:-1].replace('T', 'U')
    if len(seq):
        seqs.append(seq) 
    fp.close()
    
    return seqs, seq_targets

def get_features():
    all_weights = []
    for layer in model.layers:
       w = layer.get_weights()
       all_weights.append(w)
       
    return all_weights

def convout1_f(X):
    # The [0] is to disable the training phase flag
    return _convout1_f([0] + [X])

def get_feature(model, X_batch, index):
    inputs = [K.learning_phase()] + [model.inputs[index]]
    _convout1_f = K.function(inputs, model.layers[0].layers[index].layers[1].output)
    activations =  _convout1_f([0] + [X_batch[index]])
    
    return activations

def get_motif(model, testing, protein, y, index = 0, dir1 = 'seq_cnn/', structure  = None):
    sfilter = model.layers[0].layers[index].layers[0].get_weights()
    filter_weights_old = np.transpose(sfilter[0][:,0,:,:], (2, 1, 0)) #sfilter[0][:,0,:,:]
    print filter_weights_old.shape
    #pdb.set_trace()
    filter_weights = []
    for x in filter_weights_old:
        #normalized, scale = preprocess_data(x)
        #normalized = normalized.T
        #normalized = normalized/normalized.sum(axis=1)[:,None]
        x = x - x.mean(axis = 0)
        filter_weights.append(x)
        
    filter_weights = np.array(filter_weights)
    #pdb.set_trace()
    filter_outs = get_feature(model, testing, index)
    #pdb.set_trace()
    
    #sample_i = np.array(random.sample(xrange(testing.shape[0]), 500))
    sample_i =0

    out_dir = dir1 + protein
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    if index == 0:    
        get_motif_fig(filter_weights, filter_outs, out_dir, protein, sample_i)
    else:
        get_structure_motif_fig(filter_weights, filter_outs, out_dir, protein, y, sample_i, structure)
    
def run_network(model, total_hid, training, testing, y, validation, val_y, protein=None, structure = None):
    model.add(Dense(2, input_shape=(total_hid,)))
    model.add(Activation('softmax'))
    
    #sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
    model.compile(loss='categorical_crossentropy', optimizer='rmsprop')
    #pdb.set_trace()
    print 'model training'
    #checkpointer = ModelCheckpoint(filepath="models/" + protein + "_bestmodel.hdf5", verbose=0, save_best_only=True)
    earlystopper = EarlyStopping(monitor='val_loss', patience=5, verbose=0)

    model.fit(training, y, batch_size=50, nb_epoch=30, verbose=0, validation_data=(validation, val_y), callbacks=[earlystopper])
    
    #pdb.set_trace()
    #get_motif(model, testing, protein, y, index = 0, dir1 = 'seq_cnn1/')
    #get_motif(model, testing, protein, y, index = 1, dir1 = 'structure_cnn1/', structure = structure)

    predictions = model.predict_proba(testing)[:,1]
    return predictions, model

def run_randomforest_classifier(data, labels, test):
    clf = RandomForestClassifier(n_estimators=50)
    clf.fit(data, labels)
    #pdb.set_trace()
    pred_prob = clf.predict_proba(test)[:,1]
    return pred_prob, clf  

def run_svm_classifier(data, labels, test):
    #C_range = 10.0 ** np.arange(-1, 2)
    #param_grid = dict(C=C_range.tolist())
    clf = svm.SVC(probability =True, kernel = 'linear')
    #grid = GridSearchCV(svr, param_grid)
    clf.fit(data, labels)
    
    #clf = grid.best_estimator_
    pred_prob = clf.predict_proba(test)[:,1]
    return pred_prob, clf
    
def calculate_auc(net, hid, train, test, true_y, train_y, rf = False, validation = None, val_y = None, protein = None, structure = None):
    #print 'running network' 
    if rf:
        print 'running oli'
        #pdb.set_trace()
        predict, model = run_svm_classifier(train, train_y, test)
    else:
        predict, model = run_network(net, hid, train, test, train_y, validation, val_y, protein = protein, structure = structure)

    
    auc = roc_auc_score(true_y, predict)
    
    print "Test AUC: ", auc
    return auc, predict



def run_seq_struct_cnn_network(protein, seq = True, fw = None, oli = False, min_len = 301):
    training_data = load_data("./datasets/clip/%s/30000/training_sample_0" % protein, seq = seq, oli = oli)
    
    seq_hid = 16
    struct_hid = 16
    #pdb.set_trace()
    train_Y = training_data["Y"]
    print len(train_Y)
    #pdb.set_trace()
    training_indice, training_label, validation_indice, validation_label = split_training_validation(train_Y)
    #pdb.set_trace()
    if seq:
        cnn_train  = []
        cnn_validation = []
        seq_data = training_data["seq"][0]
        #pdb.set_trace()
        seq_train = seq_data[training_indice]
        seq_validation = seq_data[validation_indice] 
        struct_data = training_data["seq"][1]
        struct_train = struct_data[training_indice]
        struct_validation = struct_data[validation_indice] 
        cnn_train.append(seq_train)
        cnn_train.append(struct_train)
        cnn_validation.append(seq_validation)
        cnn_validation.append(struct_validation)        
        seq_net =  get_cnn_network()
        seq_data = []
            
    y, encoder = preprocess_labels(training_label)
    val_y, encoder = preprocess_labels(validation_label, encoder = encoder) 
    
    training_data.clear()
    
    rf = False

    test_data = load_data("./datasets/clip/%s/30000/test_sample_0" % protein, seq = seq, oli = oli)
    print len(test_data)
    true_y = test_data["Y"].copy()
    
    print 'predicting'    
    if seq:
        seq_test = test_data["seq"]
        structure = test_data["structure"]
        seq_auc, seq_predict = calculate_auc(seq_net, seq_hid + struct_hid, cnn_train, seq_test, true_y, y, validation = cnn_validation,
                                              val_y = val_y, protein = protein,  rf= rf, structure = structure)
        seq_train = []
        seq_test = []
         
        
        
    print str(seq_auc)
    fw.write( str(seq_auc) +'\n')

    mylabel = "\t".join(map(str, true_y))
    myprob = "\t".join(map(str, seq_predict))  
    fw.write(mylabel + '\n')
    fw.write(myprob + '\n')

def split_training_validation(classes, validation_size = 0.2, shuffle = False):
    """split sampels based on balnace classes"""
    num_samples=len(classes)
    classes=np.array(classes)
    classes_unique=np.unique(classes)
    num_classes=len(classes_unique)
    indices=np.arange(num_samples)
    #indices_folds=np.zeros([num_samples],dtype=int)
    training_indice = []
    training_label = []
    validation_indice = []
    validation_label = []
    for cl in classes_unique:
        indices_cl=indices[classes==cl]
        num_samples_cl=len(indices_cl)

        # split this class into k parts
        if shuffle:
            random.shuffle(indices_cl) # in-place shuffle
        
        # module and residual
        num_samples_each_split=int(num_samples_cl*validation_size)
        res=num_samples_cl - num_samples_each_split
        
        training_indice = training_indice + [val for val in indices_cl[num_samples_each_split:]]
        training_label = training_label + [cl] * res
        
        validation_indice = validation_indice + [val for val in indices_cl[:num_samples_each_split]]
        validation_label = validation_label + [cl]*num_samples_each_split

    training_index = np.arange(len(training_label))
    random.shuffle(training_index)
    training_indice = np.array(training_indice)[training_index]
    training_label = np.array(training_label)[training_index]
    
    validation_index = np.arange(len(validation_label))
    random.shuffle(validation_index)
    validation_indice = np.array(validation_indice)[validation_index]
    validation_label = np.array(validation_label)[validation_index]    
    
            
    return training_indice, training_label, validation_indice, validation_label        
        

def plot_roc_curve(labels, probality, legend_text, auc_tag = True):
    #fpr2, tpr2, thresholds = roc_curve(labels, pred_y)
    fpr, tpr, thresholds = roc_curve(labels, probality) #probas_[:, 1])
    roc_auc = auc(fpr, tpr)
    if auc_tag:
        rects1 = plt.plot(fpr, tpr, label=legend_text +' (AUC=%6.3f) ' %roc_auc)
    else:
        rects1 = plt.plot(fpr, tpr, label=legend_text )

def read_protein_name(filename='proteinnames'):
    protein_dict = {}
    with open(filename, 'r') as fp:
        for line in fp:
            values = line.rstrip('\r\n').split('\t')
            key_name = values[0][1:-1]
            protein_dict[key_name] = values[1]
    return protein_dict
    
def read_result_file(filename = 'result_file_seq_wohle_new'):
    results = {}
    with open(filename, 'r') as fp:
        index = 0
        #protein = '28'
        for line in fp:
            values = line.rstrip('\r\n').split('\t')
            if index % 3 == 0:
                protein = values[0].split('_')[0]
            if index % 3 != 0:
                results.setdefault(protein, []).append(values)
                
                
            index = index + 1
    
    return results

def read_individual_auc(filename = 'result_file_all_new'):
    results = {}
    with open(filename, 'r') as fp:
        index = 0
        #protein = '28'
        for line in fp:
            values = line.rstrip('\r\n').split('\t')
            pro = values[0].split('_')[0]
            results[int(pro)] = values[1:-1]    
    
    return results

def read_ideep_auc(filename = 'result_mix_auc_new'):
    results = {}
    with open(filename, 'r') as fp:
        index = 0
        #protein = '28'
        for line in fp:
            values = line.rstrip('\r\n').split('\t')
            #pdb.set_trace()
            pro = values[0].split('_')[0]
            results[int(pro)] = values[1]  
    
    return results

def plot_ideep_indi_comp():
    proteins = read_protein_name()
    ideep_resut = read_ideep_auc(filename='result_mix_auc_new')
    #pdb.set_trace()
    indi_result = read_individual_auc()
    keys = indi_result.keys()
    keys.sort()
    
    new_results = []
    names = []
    for key in keys:
        str_key = str(key)
        names.append(proteins[str_key])
        tmp = []
        for val in indi_result[key]:
            tmp.append(float(val))
        #for val in ideep_resut[key]:
        tmp.append(float(ideep_resut[key]))
        #tmp = indi_result[key] + ideep_resut[key]
        new_results.append(tmp)
    pdb.set_trace()
    new_results = map(list, zip(*new_results))
    #plot_confusion_matrix(new_results)
    plot_parameter_bar(new_results, names)
            
def plot_figure():
    protein_dict = read_protein_name()
    results = read_result_file()
    
    Figure = plt.figure(figsize=(12, 15))
    
    for key, values in results.iteritems():
        protein = protein_dict[key]
        #pdb.set_trace()
        labels = [int(float(val)) for val in values[0]]
        probability = [float(val) for val in values[1]]
        plot_roc_curve(labels, probability, protein)
    #plot_roc_curve(labels[1], probability[1], '')
    #plot_roc_curve(labels[2], probability[2], '')
    
    #title_type = 'stem cell circRNAs vs other circRNAs'
    title_type = 'ROC'
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    #plt.title(title_type)
    plt.legend(loc="lower right")
    plt.savefig('roc1.eps', format='eps') 
    #plt.show() 
  
def read_fasta_file(fasta_file):
    seq_dict = {}    
    fp = gzip.open(fasta_file, 'r')
    name = ''
    name_list = []
    for line in fp:
        line = line.rstrip()
        #distinguish header from sequence
        if line[0]=='>': #or line.startswith('>')
            #it is the header
            name = line[2:] #discarding the initial >
            name_list.append(name)
            seq_dict[name] = ''
        else:
            seq_dict[name] = seq_dict[name] + line.upper().replace('U', 'T')
    fp.close()
    
    return seq_dict, name_list

            
def get_binding_motif_fea():
    data_dir = '/home/panxy/eclipse/ideep/datasets/clip/'
    protein_dirs = os.listdir(data_dir)
    for protein in protein_dirs:
        new_dir = data_dir + protein + '/30000/'
        for beddir in os.listdir(new_dir):
            print beddir
            #pdb.set_trace()
            path = new_dir + beddir
            fas_name = os.path.join(path, 'sequences.fa')
            fw = open(fas_name, 'w')
            seq_file = os.path.join(path, 'sequences.fa.gz')
            seq_dict, name_list = read_fasta_file(seq_file)
            
            for name in name_list:
                values = name.rstrip().split(';')
                posi = values[0].split(',')
                coor = posi[0] + ':' + posi[2] + '-' + posi[3]
                fw.write('>' + coor + '\n')
                fw.write(seq_dict[name] + '\n') 
            fw.close()
            
            #pdb.set_trace()
            clistr = './get_RNA_motif_fea.sh ' + path + '>/dev/null 2>&1'
            f_cli = os.popen(clistr, 'r')
            f_cli.close()

def run_predict():
    data_dir = './datasets/clip'
    fw = open('result_file_struct_auc', 'w')
    for protein in os.listdir(data_dir):
        print protein

        fw.write(protein + '\t')

        run_seq_struct_cnn_network(protein, seq = True, fw= fw)

    fw.close()

def calculate_perofrmance(inputfile='../comp_result'):
    result = []
    with open(inputfile) as fp:
        for line in fp:
            values = line.rstrip().split('&')
            result.append([float(values[0]), float(values[1])])
    res_array = numpy.array(result)
    print np.mean(res_array, axis=0)
    print np.std(res_array, axis=0)

def read_meta_file(inputfile):
    fp =open(inputfile, 'r')
    head = True
    fw = open('download.sh', 'w')
    for line in fp:
        if head:
            head = False
            continue
        if '16 tissues mixture' in line:
            continue
        values = line.rstrip().split()
        tissue = values[-1]
        url = values[-7]
        outname = url.split('/')[-1].split('.')[0]
        new_cli = 'curl "' + url + '" -o "' + outname + '_' + tissue + '.fastq.gz"'
        fw.write(new_cli + '\n')
    fp.close()
    fw.close()
             
if __name__ == "__main__":
    run_predict()

