#!/usr/bin/env python

# Illustrates the SALIGN multiple structure/sequence alignment

import sys,os,getopt
from modeller import *

if len(sys.argv) != 5:
    #raise ValueError("Please provide a PDB accession list, one per line in capital letters with an under separating the chain (e.g. 4CMN_A)")
    print("Usage: %s [-p pdb_acc_list] [-n numb_templates]\n" % sys.argv[0])
    sys.exit(2)
def main(argv):
   infile = ''
   ntemp = ''
   try:
       opts, args = getopt.getopt(argv,"hp:n:",["pdbfile=","ntemp="])
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
      elif opt in ("-n", "--ntemp"):
         ntemp = int(arg)
         
   print('Input file is "', infile)
   
   pdb_acc = open(infile,"r").read()
   accs = pdb_acc.split("\n")
   accs_list = []
   accs_dic = {}
   for acc in accs:
       ac = tuple(acc.split("_"))
       ac_base = ac[0]
       ac_dic = {ac_base:ac}
       accs_dic.update(ac_dic)
   acc_input = list(accs_dic.values())
   
   #--- Set working dir
   p_dir = os.path.dirname(os.path.abspath(__name__))
   pdb_fpath = os.path.join(p_dir,'..')
   #pdb_accs = os.path.join(pdb_fpath,acc_input)
   
   log.verbose()
   env = environ()
   env.io.atom_files_directory = './:../atom_files/'
   
   aln = alignment(env)
   for (code, chain) in (acc_input[:ntemp]):
       mdl = model(env, file=code, model_segment=('FIRST:'+chain, 'LAST:'+chain))
       aln.append_model(mdl, atom_files=code, align_codes=code+chain)
   
   #('4cmn', 'A'), ('3qis', 'A'), ('2qv2', 'A'), ('3mtc', 'A'), 
   #        ('3n9v', 'A'), ('3qbt', 'B'), ('2kie', 'A'), ('1i9y', 'A'), 
   #        ('3nr8', 'A'), ('5okm', 'A')
   
   for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                       ((1., 0.5, 1., 1., 1., 0.), False, True),
                                       ((1., 1., 1., 1., 1., 0.), True, False)):
       aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
                  rr_file='$(LIB)/as1.sim.mat', overhang=30,
                  gap_penalties_1d=(-450, -50),
                  gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
                  dendrogram_file='fm00495.tree',
                  alignment_type='tree', # If 'progresive', the tree is not
                                         # computed and all structues will be
                                         # aligned sequentially to the first
                  feature_weights=weights, # For a multiple sequence alignment only
                                           # the first feature needs to be non-zero
                  improve_alignment=True, fit=True, write_fit=write_fit,
                  write_whole_pdb=whole, output='ALIGNMENT QUALITY')
   
   aln.write(file='fm00495.pap', alignment_format='PAP')
   aln.write(file='fm00495.ali', alignment_format='PIR')
   
   aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
              rr_file='$(LIB)/as1.sim.mat', overhang=30,
              gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
              gap_gap_score=0, gap_residue_score=0, dendrogram_file='1is3A.tree',
              alignment_type='progressive', feature_weights=[0]*6,
              improve_alignment=False, fit=False, write_fit=True,
              write_whole_pdb=False, output='QUALITY')

if __name__ == "__main__":
   main(sys.argv[1:])
