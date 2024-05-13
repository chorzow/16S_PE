###############################################
##Dmitry Sutormin, 2019##
##16S data analysis##

#1) Script parses MMseqs2 clustering results, identifies clusters with metagenome-derived sequences.
#2) Outputs clusters containing ASVs and ASVs.
###############################################


#######
#Packages to be imported.
#######
import os
import sys
import csv
from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd

#######
#Data to be used.
#######

# Path to the working directory from the pipeline settings
PWD = sys.argv[1]

#Path to the working directory.

#MMseqs2 clustering results.
MMseq_data_inpath = os.path.join(PWD, "mmseq/new_clu.tsv")
#Specify identity level used for clastering by MMseqs2.
MMseqs2_identity_level = 98
#ASV multifasta.
ASV_fasta_inpath = os.path.join(PWD, "seqtabnochim.fasta")
#ASV frequences.
ASV_freq_table_inpath = os.path.join(PWD, "seqtabnochim_n.tsv")

#Output:
Output_path = os.path.join(PWD, "OTUs_output/")


#######
#Read clustering results, prepare dictionary of clusters.
#######


def read_clusters_data(MMseq_datapath):
    with open(MMseq_datapath, 'r') as filein:
        Clusters_dict = {}
        for line in filein:
            line = line.rstrip('\r\n').split('\t')
            cluster_rep = line[0]
            cluster_dep = line[1]
            #print(cluster_rep, cluster_dep)
            if cluster_rep in Clusters_dict:
                current_cluster_ar = Clusters_dict[cluster_rep]
                current_cluster_ar.append(cluster_dep)
                Clusters_dict[cluster_rep] = current_cluster_ar
                #print('Update cluster ', Clusters_dict[cluster_rep])
            else:
                Clusters_dict[cluster_rep] = [cluster_dep]
                #print('New cluster ', Clusters_dict[cluster_rep])
    
    print('Number of clusters: ' + str(len(Clusters_dict)))
    return Clusters_dict


#######
#Analyse dictionary of clusters: clusters size distribution.
#######

def analyse_clusters_dict(Clusters_dict, Identity_level, Outpath):
    Clusters_len = []
    for cluster_rep, cluster_dep_ar in Clusters_dict.items():
        Clusters_len.append(len(cluster_dep_ar))
    
    #Plot distribution of clusters size.
    plt.hist(Clusters_len)
    plt.yscale('log')
    plt.xlabel(f'Size of clusters, {Identity_level}% identity')
    plt.ylabel('Number of clusters')
    plt.savefig(Output_path + 'MMseqs2_clusters_size_distribution.png')

    return


#######
#Read fasta with ASV sequences.
#######


def read_mfa(ASV_fasta_path):
    with open(ASV_fasta_path, 'r') as filein:
    #Keep data in dictionary.
        ASV_dict = {}
        for record in SeqIO.parse(filein, "fasta"):
            seq_16S = str(record.seq)
            id_16S = record.name
            ASV_dict[id_16S] = seq_16S
    return ASV_dict


#######
#Sum frequences of ASV comprising a cluster. Return a cumulative frequency for representative sequence.
#######

def rep_seq_cumulative_frequency(ASV_input_table_path, MMseq_Clusters_dict, ASV_sequences_dict, Outpath):
    ASV_freq_dataframe = pd.read_csv(ASV_input_table_path, sep='\t', header=0, index_col=0)
    print('Shape of the initial ASV table: ', ASV_freq_dataframe.shape)

    Rep_seq_freq_dataframe = pd.DataFrame(columns=ASV_freq_dataframe.columns)
    #Iterate clusters of sequences and sum frequences. Report cumulative frequency for representative ASV, named thereafter OTU.
    for Rep_Seq_ID, Dep_Seq_ar in MMseq_Clusters_dict.items():
        print(Dep_Seq_ar)
        cluster_freq_df = ASV_freq_dataframe.loc[Dep_Seq_ar,:]
        cluster_freq_series = cluster_freq_df.sum()
        OTU_seq_series = pd.Series([ASV_sequences_dict[Rep_Seq_ID]])
        OTU_seq_series.index = ['OTU_sequence']
        cluster_freq_series_seq = pd.concat([OTU_seq_series, cluster_freq_series])
        cluster_freq_series_seq.name = Rep_Seq_ID
        Rep_seq_freq_dataframe = Rep_seq_freq_dataframe.append(cluster_freq_series_seq)

    print('Shape of the resultant OTU table: ', Rep_seq_freq_dataframe.shape)
    Column_order = ['OTU_sequence'] + list(ASV_freq_dataframe.columns)
    Rep_seq_freq_dataframe.to_csv(path_or_buf=Outpath + "OTU_table.tsv",
                                  sep='\t', quoting=csv.QUOTE_NONNUMERIC,
                                  columns=Column_order, header=True, index=True)
    return Rep_seq_freq_dataframe


#######
#Wrapper function.
#######


def Wrapper_func_ASV_to_OTU(MMseq_datapath, Identity_level, ASV_fasta_path, ASV_input_table_path, Outpath):
    #Import mmseqs2 output.
    MMseq_Clusters_dict = read_clusters_data(MMseq_datapath)
    #Analyse distribution of cluster sizes.
    analyse_clusters_dict(MMseq_Clusters_dict, Identity_level, Outpath)
    #Import fasta file with ASV sequences.
    ASV_sequences_dict = read_mfa(ASV_fasta_path)
    #Report cumulative frequency fo OTUs.
    Rep_seq_freq_dataframe = rep_seq_cumulative_frequency(ASV_input_table_path,
                                                          MMseq_Clusters_dict,
                                                          ASV_sequences_dict,
                                                          Outpath)
    return

Wrapper_func_ASV_to_OTU(MMseq_data_inpath,
                        MMseqs2_identity_level,
                        ASV_fasta_inpath,
                        ASV_freq_table_inpath,
                        Output_path)
