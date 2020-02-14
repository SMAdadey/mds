#!/usr/bin/env python

import sys,os,getopt
from modeller import *
from modeller.automodel import *

if len(sys.argv) != 9:
    #raise ValueError("Please provide a PDB accession list, one per line in capital letters with an under separating the chain (e.g. 4CMN_A)")
    print("Usage: %s [-p pdb_acc_list] [-s seq_pir_file] [-n numb_templates] [-m num_models]\n" % sys.argv[0])
    usg = """Build Models
        
        -p,--pdbfile   <str>      :PDB accessions list (e.g. blast_result.acc)
        -s,--pirfile   <str>      :Your sequence file in PIR format (e.g. gjb4.ali)
        -n,--ntemp     <int>      :Number of templates to use
        -m,--nmodel    <int>      :Number of models to build
    """
    print(usg)
    sys.exit(2)
def main(argv):
   infile = ''
   ntemp = ''
   pirfile = ''
   nmod = ''
   try:
       opts, args = getopt.getopt(argv,"hp:s:n:m:",["pdbfile=","pirfile=","ntemp=","nmodel="])
   except getopt.GetoptError:
      print('%s -p <pdbfile_list> -s <seq_pir_file> -n <num_templates> -m <num_models>' % sys.argv[0])
      sys.exit(2)
   for opt, arg in opts:
      if opt in ('-h','--help'):
         print('%s -p <pdbfile_list> -s <seq_pir_file> -n <num_templates> -m <num_models>' % sys.argv[0])
         sys.exit()
      elif opt in (" ", ""):
         print('%s -p <pdbfile_list> -s <seq_pir_file> -n <num_templates> -m <num_models>' % sys.argv[0])
         sys.exit()
      elif opt in ("-p", "--pdbfile"):
         infile = arg
      elif opt in ("-s", "--pirfile"):
         pirfile = arg
      elif opt in ("-n", "--ntemp"):
         ntemp = int(arg)
      elif opt in ("-m", "--nmodel"):
         nmod = int(arg)

   print('Input file is "', infile)
   pir_base = pirfile.replace(".ali","")

   pdb_acc = open(infile,"r").read()
   accs = pdb_acc.split("\n")
   accs_list = []
   accs_dic = {}
   for acc in accs:
       ac = acc.replace("_","")
       acc_b = acc.split("_")
       ac_base = acc_b[0]
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
   a.ending_model = nmod
   a.make()

if __name__ == "__main__":
   main(sys.argv[1:])
