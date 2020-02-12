from modeller import *
import sys,os,getopt

if len(sys.argv) != 3:
    #raise ValueError("Please provide a PDB accession list, one per line in capital letters with an under separating the chain (e.g. 4CMN_A)")
    print("Usage: %s [-p pir_seq_file]\n" % sys.argv[0])
    sys.exit(2)
def main(argv):
   infile = ''
   try:
       opts, args = getopt.getopt(argv,"hp:",["pirfile="])
   except getopt.GetoptError:
      print('align2d.py -p <pirfile=>')
      sys.exit(2)
   for opt, arg in opts:
      if opt in ('-h','--help'):
         print('align2d.py -p <pirfile>')
         sys.exit()
      elif opt in (" ", ""):
         print('align2d.py -p <pirfile>')
         sys.exit()
      elif opt in ("-p", "--pdbfile"):
         infile = arg

   print('Input file is "', infile)
   pir_base = infile.replace(".ali","")

   log.verbose()
   env = environ()
   
   env.libs.topology.read(file='$(LIB)/top_heav.lib')
   
   # Read aligned structure(s):
   aln = alignment(env)
   aln.append(file='fm00495.ali', align_codes='all')
   aln_block = len(aln)
   
   # Read aligned sequence(s):
   aln.append(file='%s.ali' % pir_base.lower(), align_codes='%s' % pir_base.upper())
   
   # Structure sensitive variable gap penalty sequence-sequence alignment:
   aln.salign(output='', max_gap_length=20,
              gap_function=True,   # to use structure-dependent gap penalty
              alignment_type='PAIRWISE', align_block=aln_block,
              feature_weights=(1., 0., 0., 0., 0., 0.), overhang=0,
              gap_penalties_1d=(-450, 0),
              gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
              similarity_flag=True)
   
   aln.write(file='%s-mult.ali' % pir_base.lower(), alignment_format='PIR')
   aln.write(file='%s-mult.pap' % pir_base.lower(), alignment_format='PAP')
if __name__ == "__main__":
    main(sys.argv[1:])
