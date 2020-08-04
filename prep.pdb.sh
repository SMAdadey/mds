#!/usr/bin/env bash

#Prepare PDB templates for Modelling: Extract the chains required

if [ $# != 1 ]; then
   echo "Usage: prep.pdb.sh [file.acc]"
   echo "NB: Make sure PDB files listed in file.acc are present in the same directory as file.acc"
else
   acc="$1"
   for i in $(cat $acc); do 
     if [[ -f "${i/_*/.pdb}" ]]; then
       mkdir -p raw_pdbs; 
       mv ${i/_*/.pdb} raw_pdbs/;
     fi
       echo "-c -d \'load raw_pdbs/${i/_*/.pdb}; indicate chain ${i/*_/}; save ${i/_*/}.pdb\'"; 
   done | \
        parallel echo pymol {} | xargs -I input -P11 sh -c "input"
fi
