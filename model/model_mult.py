#!/usr/bin/env python

import sys,os,getopt
from modeller import *
from modeller.automodel import *

if len(sys.argv) != 7:
    #raise ValueError("Please provide a PDB accession list, one per line in capital letters with an under separating the chain (e.g. 4CMN_A)")
    print("Usage: %s [-p pdb_acc_list] [-s pir_seq_file] [-n numb_templates]\n" % sys.argv[0])
    sys.exit(2)
def main(argv):
   infile = ''
   pirfile = ''
   ntemp = ''
   try:
       opts, args = getopt.getopt(argv,"hp:s:n:",["pdbfile=","pirfile=","ntemp="])
   except getopt.GetoptError:
      print('salign.py -p <pdbfile_list> -n <num_templates>')
      sys.exit(2)
   for opt, arg in opts:
      if opt in ('-h','--help'):
         print('salign.py -p <pdbfile_list> -n <num_templates>')
         sys.exit()
      elif opt in (" ", ""):
         print('salign.py -p <pdbfile_list> -n <num_templates>')
         sys.exit()
      elif opt in ("-p", "--pdbfile"):
         infile = arg
      elif opt in ("-s", "--pirfile"):
         pirfile = arg
      elif opt in ("-n", "--ntemp"):
         ntemp = int(arg)

   print('Input file is "', infile)
   pir_base = pirfile.replace(".ali","")

   pdb_acc = open(infile,"r").read()
   accs = pdb_acc.split("\n")
   accs_list = []
   accs_dic = {}
   for acc in accs:
       ac = acc.replace("_","")
       ac_base = ac[0]
       ac_dic = {ac_base:ac}
       accs_dic.update(ac_dic)
       accs_list.append(ac)
   acc_input = list(accs_list)
   print(accs_list)
   env = environ()
   a = automodel(env, alnfile='%s-mult.ali' % pir_base.lower(),
                 knowns=(acc_input[:ntemp]), sequence='%s' % pir_base.upper(), 
                 assess_methods=(assess.DOPE,
                                 #soap_protein_od.Scorer(),
                                 assess.GA341))
   
   #('4cmnA', '3qisA', '2qv2A', 
   #              '3mtcA', '3n9vA', '3qbtB', '2kieA', 
   #              '1i9yA', '3nr8A', '5okmA')
   
   a.starting_model = 1
   a.ending_model = 50
   a.make()

if __name__ == "__main__":
   main(sys.argv[1:])
