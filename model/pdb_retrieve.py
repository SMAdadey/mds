#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: Tuesday Feb 11, 2020

@author: Kevin Esoh
"""

import sys,os,os.path, getopt
import urllib.request as ur
from Bio.Blast import NCBIWWW as ncbiweb
from Bio.Blast import NCBIXML as ncbiparser

if len(sys.argv) != 7:
    #raise ValueError("Please provide a PDB accession list, one per line in capital letters with an under separating the chain (e.g. 4CMN_A)")
    print("Usage: %s [-s seq_file] [-o out_file] [-n num_seq]\n" % sys.argv[0])
    sys.exit(2)
def main(argv):
   #global infile 
   #global outfile
   infile = ''
   outfile = ''
   nseq = ''
   try:
       opts, args = getopt.getopt(argv,"hs:o:n:",["seqfile=","outfile=","nseq="])
   except getopt.GetoptError:
      print('pdb_retrieve.py -s <seqlist> -o <out_file>')
      sys.exit(2)
   for opt, arg in opts:
      if opt in ('-h','--help'):
         print('pdb_retrieve.py -s <seqlist> -o <out_list>')
         sys.exit()
      elif opt in (" ", ""):
         print('pdb_retrieve.py -s <seqlist> -o <out_list>')
         sys.exit()
      elif opt in ("-s", "--seqfile"):
         infile = arg
      elif opt in ("-o", "--outfile"):
         outfile = arg
      elif opt in ("-n", "--nseq"):
         nseq = int(arg)


   print('Input file is "', infile)

   #--- load query sequence
   # s = 'ocrl_mt_prot.fasta'
   seq = open(infile,"r").read()
    
   # #--- open file to write out blast xml results
   # outfile = "blast_result.xml"
   blast_result = open(outfile, "w+")
    
   #--- blast search
   print("BLAST search running...")
   blast_handle = ncbiweb.qblast(program="blastp", database="pdbaa", sequence=seq, service="psi", hitlist_size=nseq) # expect=10
   blast_hits = blast_handle.read()
   blast_result.write(blast_hits)
   blast_result.close()
   print("BLAST completed successful!")
   
   #--- parse blast results
   blast_hits = ncbiparser.parse(open(outfile))    # blast_hits is a generator containing records
   
   hit_acc = []    # get accessions
   pdb_input = []   # make PDB input from accessions
   acclib = {}
   pdb_acclib = {} # make a dict of PDB accessions and their chains
   acclist = []
   acc_file = open('blast_result.acc','w+')
   
   #--------------------- NEW ------- NEW ----------- NEW -------------------
   hits_acc = []
   hits_tup = ()
   hits_dic = {}
   acclist = []
   for record in blast_hits:
       for hit in record.alignments:
           acclist.append(hit.accession)
           acc = hit.accession.split("_")
           hits_tup = tuple(acc)
           acc_base = acc[0]
           accdic = {acc_base:hits_tup}
           hits_dic.update(accdic)
           hits_acc.append("%s.pdb" % acc_base)
           #print("PDB INPUT: %s.pdb" % acc_base.lower())
   acc_file.write('\n'.join(acclist))
   acc_file.close()
   print("<<<< >>>>")
   print("PDB ACCs: %s\n" % acclist)
   #print(pdb_acclib['4cmn'])
   #------------------ NEW ----------- NEW ------------- NEW ---------------
   
   # for record in blast_hits:
   #    for hit in record.alignments:
   #        acclist.append(hit.accession)
   #        acc = hit.accession.split("_")
   #        hit_acc.append(acc)
   #        acc_base = acc[0].lower()
   #        acclib = {acc_base:acc}
   #        pdb_acclib.update(acclib)
   #        pdb_input.append("%s.pdb" % acc_base)
   # acc_file.write('\n'.join(acclist))
   # acc_file.close()
   # print("<<<< >>>>")
   # print("PDB files: %s\n" % pdb_input)
   # #print(pdb_acclib['4cmn'])
   
   #--- download PDB coordinate files
   print("Fetching pdb files...")
   for acc in hits_acc:
       #ac = acc.split(".")
       if os.path.exists(acc):
           print("%s already exists! Skipping..." % acc)
           continue
       else :
           ur.urlretrieve("http://files.rcsb.org/download/%s" % acc.lower(), acc)    
   print("")
   print("Done fetching pdb files!")

if __name__ == "__main__":
   main(sys.argv[1:])
