from eden.converter.fasta import sequence_to_eden
from eden.modifier.rna.annotate_rna_structure import annotate_single
import subprocess as sp
import os
import pdb
import gc

def run_rnashape(sequence):
    #
    cmd = 'echo "%s" | ./RNAshapes -t %d -c %d -# %d' % (sequence,5, 10, 1)
    #out = sp.check_output(cmd, shell=True)
    #print out
    text = []
    fp = os.popen(cmd, 'r')
    for line in fp:
        text.append(line)
    fp.close()
    #text = out.strip().split('\n')
    #out.close()
    seq_info = text[0]
    if 'configured to print' in text[-1]:
        struct_text = text[-2]
    else:
        struct_text = text[1]
    # shape:
    structur = struct_text.split()[1]
    # extract the shape bracket notation
    #shape_list += [line.split()[2] for line in struct_text]
    #encoee strucyrte
    graph = sequence_to_eden([("ID", sequence)]).next()
    graph.graph['structure']=structur
    annotate_single(graph)
    encode_struct = ''.join([ x["entity_short"].upper() for x in graph.node.values() ])
    gc.collect()
    return encode_struct
    #pdb.set_trace()

def read_structure(seq_file):
    seq_list = []
    structure_list = []
    fw = open(seq_file + '.structure', 'w')
    seq = ''
    with open(seq_file, 'r') as fp:
        for line in fp:
            if line[0] == '>':
                name = line
                if len(seq):
                    fw.write(old_name)
                    seq = seq.replace('U', 'T')
                    struc_en = run_rnashape(seq)
                    fw.write(struc_en + '\n')
                old_name = name              
                seq = ''
            else:
                seq = seq + line[:-1]
        if len(seq):
            fw.write(old_name)
            seq = seq.replace('U', 'T')
            struc_en = run_rnashape(seq)
            fw.write(struc_en + '\n') 
    fw.close()

def run_predict_graphprot_data():
    data_dir = '/home/panxy/eclipse/ideep/GraphProt_CLIP_sequences/'
    #fw = open('result_file_struct_graphprot', 'w')
    for protein_file in os.listdir(data_dir):

        protein = protein_file.split('.')[0]
        print protein
        read_structure(data_dir + protein_file)
        #fw.write(protein + '\t')
        #model = merge_seperate_network_with_multiple_features(protein, kmer=False, rg=True, clip=True, rna=True, go=False, motif = True, seq = True, fw = fw)
        #run_individual_network(protein, kmer=False, rg=False, clip=False, rna=False, motif = False, seq = True, fw = fw)
        #run_seq_struct_cnn_network(protein, seq = True, fw= fw, graph_prot = True)
        #run_individual_network(protein, kmer=False, rg=False, clip=False, rna=False, go=False, motif = False, seq = True, oli = False, fw = fw)
    #fw.close()

if __name__ == "__main__":  
    #run_predict_graphprot_data()  
    sequence = 'TGGAAACATTCCTCAGGTGGTTCATCCAAGGCCCTTTCCACTCTTTCAGCTCACAGCACAGTGGTCCTTTTGTTCTTTGGTCCACCCATGTTTGTGTATAC'
    encode_struct = run_rnashape(sequence)
    print encode_struct
