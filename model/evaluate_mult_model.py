#!/usr/bin/env python

import os,sys,getopt
from modeller import *
from modeller.scripts import complete_pdb

if len(sys.argv) != 7:
    print("\nUsage: %s [-p best_model_pdbfile] [-s pir_seq_file] [-t best_pdb_template]\n" % sys.argv[0])
    usg = """Evaluate Models
        
        -p,--pdbfile <str>      :Best pdb model from the modeling step (e.g. OCRL.B99990037.pdb)
        -s,--pirfile <str>      :Subject sequence file in PIR format (e.g. gjb4.ali)
        -t,--temp    <str>      :Best template pdb file (e.g. 4CMN.pdb)
    """
    print(usg)
    sys.exit(2)
def main(argv):
   infile = ''
   pirfile = ''
   ntemp = ''
   try:
       opts, args = getopt.getopt(argv,"hp:s:t:",["pdbfile=","pirfile=","temp="])
   except getopt.GetoptError:
      print('%s -p <best_model_pdbfile> -s <pir_seq_file> -t <best_template_file>' % sys.argv[0])
      sys.exit(2)
   for opt, arg in opts:
      if opt in ('-h','--help'):
         print('%s -p <best_model_pdbfile> -s <pir_seq_file> -t <best_template_file>' % sys.argv[0])
         sys.exit()
      elif opt in (" ", ""):
         print('%s -p <best_model_pdbfile> -s <pir_seq_file> -t <best_template_file>' % sys.argv[0])
         sys.exit()
      elif opt in ("-p", "--pdbfile"):
         infile = arg
      elif opt in ("-s", "--pirfile"):
         pirfile = arg
      elif opt in ("-t", "--temp"):
         btemp = arg

   print('Input file is "', infile)
   pir_base = pirfile.replace(".ali","")
   btemp_profile = btemp.replace("pdb","profile")

   log.verbose()    # request verbose output
   env = environ()
   env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
   env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters
   
   # directories for input atom files
   env.io.atom_files_directory = './:../atom_files'
   
   # read model file
   mdl = complete_pdb(env, infile)
   tpl = complete_pdb(env, btemp)
   
   s = selection(tpl)
   s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=btemp_profile,
                 normalize_profile=True, smoothing_window=15)
   
   s = selection(mdl)
   s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='%s-mult.profile' % pir_base,
                 normalize_profile=True, smoothing_window=15)

if __name__ == "__main__":
    main(sys.argv[1:])
