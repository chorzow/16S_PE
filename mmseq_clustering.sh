#!/bin/bash

echo "You can set identity and minimum coverage as second and third arguments respectively."


DIR=$1
DATAFILE="${DIR}/seqtabnochim.fasta"
DATASETPREFIX='new'
#if [[ -z $1 ]]; then
# echo "No argument has been provided :("
# echo "Please provide ouput files prefix as an argument."
# exit 1
#else
# DATASETPREFIX=$1
#fi

CPWD=$PWD

TMP_DIR="${DIR}/tmp"
MMSEQSOUTDIR="${DIR}/mmseq"

IDENTITY=0.98
if [ $2 ]; then
 IDENTITY=$2
fi
echo "Clustering with identiy: "$IDENTITY

COVERAGE=0.8
if [ $3 ]; then
 COVERAGE=$3
fi
echo "Minimum coverage: "$COVERAGE


# creating temporary files directory
if [ ! -d "$TMP_DIR" ]; then
	mkdir $TMP_DIR
fi
# creating output files directory
if [ ! -d "$MMSEQSOUTDIR" ]; then
	mkdir $MMSEQSOUTDIR
fi

cd $MMSEQSOUTDIR

# Converting FASTA file to MMSEQS database

mmseqs createdb $DATAFILE "$DATASETPREFIX"_DB

# Then execute the clustering:

mmseqs cluster "$DATASETPREFIX"_DB "$DATASETPREFIX"_clu $TMP_DIR --min-seq-id $IDENTITY -c $COVERAGE --single-step-clustering false --cluster-mode 0 -v 1 --split-memory-limit 1000000

# Please ensure that in case of large input databases the temporary direcotry provides enough free space. For disk space requirements, see the user guide.

# To generate a TSV-style formatted output file from the ffindex output file, type:
echo 'Generate DB_clu.tsv...'

mmseqs createtsv "$DATASETPREFIX"_DB "$DATASETPREFIX"_DB "$DATASETPREFIX"_clu "$MMSEQSOUTDIR"/"$DATASETPREFIX"_clu.tsv

# To generate a FASTA-style formatted output file from the ffindex output file, type:

mmseqs createseqfiledb "$DATASETPREFIX"_DB "$DATASETPREFIX"_clu "$DATASETPREFIX"_clu_seq
mmseqs result2flat "$DATASETPREFIX"_DB "$DATASETPREFIX"_DB "$DATASETPREFIX"_clu_seq "$MMSEQSOUTDIR"/"$DATASETPREFIX"_clu_seq.fasta


# To extract the representative sequences from the clustering result call:

mmseqs result2repseq "$DATASETPREFIX"_DB "$DATASETPREFIX"_clu "$DATASETPREFIX"_DB_clu_rep
mmseqs result2flat "$DATASETPREFIX"_DB "$DATASETPREFIX"_DB "$DATASETPREFIX"_DB_clu_rep "$MMSEQSOUTDIR"/"$DATASETPREFIX"_DB_clu_rep.fasta --use-fasta-header

echo "Min coverage: "$COVERAGE
cd $CPWD
rm -r $TMP_DIR

echo "Doubled cluster count"
wc -l "$MMSEQSOUTDIR"/"$DATASETPREFIX"_DB_clu_rep.fasta
