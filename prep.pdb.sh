#!/usr/bin/env bash

for i in $(cat covid19.acc); do 
  if [[ -f "${i/_*/.pdb}" ]]; then
    mkdir -p raw_pdbs; 
    mv ${i/_*/.pdb} raw_pdbs/;
  fi
    echo "-c -d \'load raw_pdbs/${i/_*/.pdb}; indicate chain ${i/*_/}; save ${i/_*/}.pdb\'"; 
done | \
     parallel echo pymol {} | xargs -I input -P11 sh -c "input"
