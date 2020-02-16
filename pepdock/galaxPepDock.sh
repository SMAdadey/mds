#!/usr/bin/env bash

function usage() {
        printf "Usage: %s [-j job_type] [-p pdb_file] [-n project_name]\n" $(basename $0)
        echo """
                Galaxy Protein-Peptide Docking

                -j,--jtype     [clus/inte]     :Submit job on cluster [clus] or run interatively [inte]
                -p,--pdb       <str>           :PDB protein stucture
                -s,--seq       <str>           :Peptide sequnce to dock in FASTA format
                -n,--name      <str>           :Name of project
                -h,--help      		       :Print this message
        """
}

if [[ $# != 8 ]]; then
   usage

else
   prog=`getopt -o "hj:p:s:n:" --long "help,jtype:,pdb:,seq:,name:" -- "$@"`
   
   if [ $? != 0 ]; then 
      echo "An error occurred! Terminating..." 1>&2; 
      exit 1; 
   fi
   
   eval set -- "$prog"
   
   while true; do
       case "$1" in
          -j|--jtype) jtp="$2"; shift 2;;
          -p|--pdb) inpdb="$2"; shift 2;;
          -s|--seq) seq="$2"; shift 2;;
          -n|--name) pname="$2"; shift 2;;
          -h|--help) shift; usage; 1>&2; exit 1 ;;
          --) shift; break ;;
          *) shift; usage; 1>&2; exit 1 ;;
       esac
   done
   
   echo Job type: $jtp
   echo "Project name: $pname"
   
   if [[ $jtp == "inte" ]]; then
      cd /mnt/lustre/groups/CBBI1243/KEVIN/git/mds/pepdock/
      export GALAXY_HOME=/mnt/lustre/users/kesoh/bioTools/GalaxyPepDock 
      echo "==========================================="
      echo "PDB file: $inpdb"
      echo "CMD: GalaxyPepDock -t $pname -p $inpdb -s $seq"
      echo "==========================================="
      sleep 1
      echo "Running..."
      GalaxyPepDock.ubuntu1604 -p $inpdb -t $pname -s $seq
   
   elif [[ $jtp == "clus" ]]; then
      echo """#!/usr/bin/env bash
#PBS -l select=1:ncpus=24
#PBS -l walltime=48:00:00
#PBS -q smp
#PBS -P CBBI1243
#PBS -o /mnt/lustre/groups/CBBI1243/KEVIN/git/mds/pepdock/stdout.txt
#PBS -e /mnt/lustre/groups/CBBI1243/KEVIN/git/mds/pepdock/stderr.txt
#PBS -N GalRefine
#PBS -M kevin.esoh@students.jkuat.ac.ke
#PBS -m b

cd /mnt/lustre/groups/CBBI1243/KEVIN/git/mds/pepdock/

set GALAXY_HOME=/mnt/lustre/users/kesoh/bioTools/GalaxyPepDock
GalaxyPepDock.ubuntu1604 -p $inpdb -t $pname -s $seq

""" > ${pname}.qsub

   echo """
       Refinement job created!
       Please submit using '${pname}.qsub'
   """
   else
      echo "Invalid option!" 1>&2;
      exit 1;
   fi
fi
