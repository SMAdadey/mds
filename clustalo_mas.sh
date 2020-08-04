#!/usr/bin/env bash

if [ $# != 3 ]; then
	echo "Usage: ./clustalo_msa.sh [in_multi_fasta] [out_format] [wrap]"
else
	in=$1; out=${in/.fa*/.aln}; of=$2; wr=$3
	sed 's/ /_/g' $in > tmp.fasta
	clustalo -i tmp.fasta -o $out --wrap=$wr --force --outfmt=$of
	if [ -e "tmp.fasta" ]; then rm tmp.fasta; fi
fi
