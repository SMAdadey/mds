#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created: Tuesday Feb 11, 2020

@author: Kevin Esoh
"""

import sys,os,os.path
import urllib.request as ur
from Bio.Blast import NCBIWWW as ncbiweb
from Bio.Blast import NCBIXML as ncbiparser

# #--- load query sequence
# s = 'ocrl_mt_prot.fasta'
# seq = open(s,"r").read()
# 
# #--- open file to write out blast xml results
outfile = "blast_result.xml"
# blast_result = open(outfile, "w+")
# 
# #--- blast search
# blast_handle = ncbiweb.qblast(program="blastp", database="pdbaa", sequence=seq, service="psi", hitlist_size=20) # expect=10
# blast_hits = blast_handle.read()
# blast_result.write(blast_hits)
# blast_result.close()

#--- parse blast results
blast_hits = ncbiparser.parse(open(outfile))    # blast_hits is a generator containing records

hit_acc = []    # get accessions
pdb_input = []   # make PDB input from accessions
acclib = {}
pdb_acclib = {} # make a library of PDB accessions and their chains

for record in blast_hits:
   for hit in record.alignments:
       acc = hit.accession.split("_")
       hit_acc.append(acc)
       acc_base = acc[0].lower()
       acclib = {acc_base:acc}
       pdb_acclib.update(acclib)
       pdb_input.append("%s.pdb" % acc_base)
print("<<<< >>>>")
print("PDB files: %s\n" % pdb_input)
#print(pdb_acclib['4cmn'])

#--- download PDB coordinate files
print("Fetching pdb files...")
for acc in pdb_input:
    if os.path.exists(acc):
        print("%s already exists! Skipping..." % acc)
        continue
    else :
        ur.urlretrieve("http://files.rcsb.org/download/%s" % acc, acc)    
print("")
print("Done fetching pdb files!")
