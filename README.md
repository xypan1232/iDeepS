# iDeepS
iDeepS aims to discover Sequence-Structure Motifs, it trains deep learning models to infer binding sequence and structure motifs from sequences simultaneously.
We first encode the sequence and secondary structure into one-hot encoding, which are further fed into CNNs to learn abstract motif features. 
then we use bidirectional LSTM to capture the long range dependencies between binding sequence and structure motifs identified by CNNs.
Finally the learned abstract features are fed into classification layer to predict RBP binding sites on RNAs.
We comprehensively evaluate iDeepS on verified RBP binding sites derived from large-scale representative CLIP-seq datasets.


# Dependency <br>
python 2.7 <br>
NumPy v1.6.1 <br>
Scipy v0.9 <br>
<a href=https://github.com/fchollet/keras/>keras v1.1.2 library</a> and its backend is theano 0.9.0 <br>
<a href=https://github.com/scikit-learn/scikit-learn>sklearn v0.17</a> <br>
<a href=https://github.com/fabriziocosta/EDeN>EDeN</a> NOTE: use our uploaded EDeN.zip, decompress it and install it locally. The code structure of latest EDeN is completely changed and it does not work with our code.<br> 
<a href=https://bibiserv.cebitec.uni-bielefeld.de/download/tools/rnashapes.html>RNAshapes</a> . If you cannot download it, please use our uploaded RNAshapes in this github repository. <br>

# Content <br>
./datasets: the training and testing dataset with sequence and label indicating it is binding sites or not<br>
./motif/seq_cnn: detected binding sequence motifs from iDeepS, and we also report the matched known motifs uusing TOMTOM and motif enrichment analysis using AME in MEME Suite. <br>
./motif/structure_cnn: detected binding structure motifs from iDeepS. It also reports the motif enrichment analysis using AME in MEME Suite<br>
./ideeps.py: the python code, it can be ran to reproduce our results. <br>


# Usage

 python ideeps.py [-h] [--data_file <data_file>] [--train TRAIN] <br>
                [--model_dir MODEL_DIR] [--predict PREDICT] <br>
                [--out_file OUT_FILE] [--batch_size BATCH_SIZE] <nr>
                [--n_epochs N_EPOCHS]  [--motif MOTIF]   [--motif_dir MOTIF_DIR] <br> <br>
where the input training file should be sequences.fa.gz with label info in each head per sequence.<br>

# Use example
<b>1.</b> Train the model using your data (currently only support fix-length sequences): <br>
python ideeps.py --train=True --data_file=datasets/clip/10_PARCLIP_ELAVL1A_hg19/30000/training_sample_0/sequences.fa.gz --model_dir=models
<br> <br>
--model_dir: the dir used to save the trained model, which is used for prediction step. <br>
 --data_file: the <b>training</b> sequence file sequences.fa.gz with label informaiton in the head. <br>
 <br>
<b>NOTICE:</b> When you run iDeepS, please make sure there are no empty structure.gz in corresponding dir. otherwise it will give the error: <br>
  File "ideeps.py", line 163, in read_structure <br>
    structure = struct_dict[old_name[:-1]] <br>
KeyError: '> chr1,+,44951749,44951849; class:0' <br>
<br>
<b>2.</b> Predict the binding probability for your sequences (you need use the same dir for saved models in training step): <br>
 python ideeps.py --predict=True --data_file=datasets/clip/10_PARCLIP_ELAVL1A_hg19/30000/test_sample_0/sequences.fa.gz --model_dir=models --out_file=YOUR_OUTFILE
<br> <br>
--model_dir: The saved dir for models in training step. <br>
--data_file: configure your <b>testing</b> sequence file sequences.fa.gz.

Need install WebLogo (http://weblogo.berkeley.edu/) and TOMTOM in MEME Suite(http://meme-suite.org/doc/download.html?man_type=web) to search identifyed motifs <br>
<b>3.</b> Identify the binding sequence-structure motifs (you need use the same dir for saved models in training step): <br>
 python ideeps.py --motif=True --data_file=datasets/clip/10_PARCLIP_ELAVL1A_hg19/30000/test_sample_0/sequences.fa.gz --model_dir=models --motif_dir=YOUR_MOTIF_DIR
<br> <br>
--model_dir: The saved dir for models in training step, you must specify the trained model. <br>
--data_file: configure your sequence file sequences.fa.gz to identify binding sequence-structure motifs.

<br><br>

# NOTE
if there is memory issue for RNA strucutre prediciton, mayeb use rnashape_structure_without_memory_issue.py instead of rnashape_structure.py. We replace sp.checkout using os.open, which maybe slower, but will not cause memory issue.

# Update
We also update iDeepS to iDeepS2, which can handle vairbale lengths and the sequences and structures are encoded into one-hot encoding vector. You can dowbload iDeepS2 from https://github.com/xypan1232/iDeepS2.

# Reference
Xiaoyong Pan^, Peter Rijnbeek, Junchi Yan, Hong-Bin Shen^. <a href = "https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4889-1">Prediction of RNA-protein sequence and structure binding preferences using deep convolutional and recurrent neural networks</a>. BMC Genomics, 2018, 19:511. <br><br>
Contact: Xiaoyong Pan (xypan172436atgmail.com)



