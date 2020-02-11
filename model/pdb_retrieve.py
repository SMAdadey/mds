#!/usr/bin/env python

import sys,os
import urllib.request as ur
from Bio.Blast import NCBIWWW as ncbiweb
from Bio.Blast import NCBIXML as ncbiparser

#--- load query sequence
s = 'ocrl_mt_prot.fasta'
seq = open(s,"r").read()

#--- open file to write out blast xml results
outfile = "blast_result.xml"
blast_result = open(outfile, "w+")

#--- blast search
blast_handle = ncbiweb.qblast(program="blastp", database="pdbaa", sequence=seq, service="psi", hitlist_size=10, expect=1)
blast_hits = blast_handle.read()
blast_result.write(blast_hits)
blast_result.close()

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
print("PDB INPUT: %s\n" % pdb_input)
#print(pdb_acclib['4cmn'])

#--- download PDB coordinate files
for acc in pdb_input:
    print("Fetching pdb files...")    
    ur.urlretrieve("http://files.rcsb.org/download/%s" % acc, acc)
    print("Done fetching pdb files!")
    print("")


