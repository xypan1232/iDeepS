# iDeepS
iDeepS aims to discover Sequence-Structure Motifs, it trains deep learning models to infer binding sequence and structure motifs from sequences simultaneously.
We first encode the sequence and secondary structure into one-hot encoding, which are further fed into CNNs to learn abstract motif features. 
then we use bidirectional LSTM to capture the long range dependencies between binding sequence and structure motifs identified by CNNs.
Finally the learned abstract features are fed into classification layer to predict RBP binding sites on RNAs.
We comprehensively evaluate DSSM on verified RBP binding sites derived from large-scale representative CLIP-seq datasets.

<br>
Dependency <br>
<a href=https://github.com/fchollet/keras/>keras 1.1.2 library</a> <br>
<a href=https://github.com/scikit-learn/scikit-learn>sklearn</a> <br>
<a href=https://github.com/fabriziocosta/EDeN>EDeN</a> <br>
<a href=https://bibiserv.cebitec.uni-bielefeld.de/download/tools/rnashapes.html>RNAshapes</a> <br>

Content <br>
./datasets: the training and testing dataset with sequence and label indicating it is binding sites or not<br>
./motif/seq_cnn: detected binding sequence motifs from iDeepS. <br>
./motif/structure_cnn: detected binding structure motifs from iDeepS. <br>
./ideeps.py: the python code, it can be ran to reproduce our results. <br>

Contact: Xiaoyong Pan (xypan172436atgmail.com)
