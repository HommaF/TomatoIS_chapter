#!/bin/bash

id_file=$1'.txt'
outfile=$1'_ms_hits.tsv'
og_outfile=$1'_out.tsv'
tryptic_outfile=$1'_tryptic_peptides.tsv'
ms_file=$2
ids=$3

head -1 $ids | cut -f1,2,3,34,35,36,39,40,43,64-74,86-93 > $outfile
cat $ms_file | xargs -i -n 1 grep -f $id_file {} | cut -f1,2,3,34,35,36,39,40,43,64-74,86-93 >> $outfile

cat $id_file | xargs -i -n 1 seqkit grep --pattern {} ~/Desktop/coding/proteome/*fasta | prot2pept > tmp_$outfile

echo $tryptic_outfile

./tryptic_detection.py $outfile tmp_$outfile $tryptic_outfile  > $og_outfile

rm tmp*
