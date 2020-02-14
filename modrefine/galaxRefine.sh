#!/usr/bin/env bash

function usage() {
        printf "Usage: %s [-j job_type] [-p pdb_file] [-n project_name]\n" $(basename $0)
        echo """
                Galaxy Model Refinement

                -j,--jtype     [clus/inte]     :Submit job on cluster [clus] or run interatively [inte]
                -p,--pdb       <str>           :PDB file to refine
                -n,--name      <str>           :Name of project
                -h,--help      		       :Print this message
        """
}

if [[ $# != 6 ]]; then
   usage

else
   prog=`getopt -o "hj:p:n:" --long "help,jtype:,pdb:,name:" -- "$@"`
   
   if [ $? != 0 ]; then 
      echo "An error occurred! Terminating..." 1>&2; 
      exit 1; 
   fi
   
   eval set -- "$prog"
   
   while true; do
       case "$1" in
          -j|--jtype) jtp="$2"; shift 2;;
          -p|--pdb) inpdb="$2"; shift 2;;
          -n|--name) pname="$2"; shift 2;;
          -h|--help) shift; usage; 1>&2; exit 1 ;;
          --) shift; break ;;
          *) shift; usage; 1>&2; exit 1 ;;
       esac
   done
   
   echo Job type: $jtp
   echo "Project name: $pname"
   
   if [[ $jtp == "inte" ]]; then
      cd /mnt/lustre/groups/CBBI1243/KEVIN/git/mds/modrefine/
      
      echo "==========================================="
      echo "PDB file: $inpdb"
      echo "CMD: GalaxyRefine -p $inpdb -t $pname -o 10"
      echo "==========================================="
      sleep 1
      echo "Running..."
      GalaxyRefine -p $inpdb -t $pname -o 10
   
   elif [[ $jtp == "clus" ]]; then
      echo """#!/usr/bin/env bash
#PBS -l select=1:ncpus=24
#PBS -l walltime=48:00:00
#PBS -q smp
#PBS -P CBBI1243
#PBS -o /mnt/lustre/groups/CBBI1243/KEVIN/git/mds/modrefine/stdout.txt
#PBS -e /mnt/lustre/groups/CBBI1243/KEVIN/git/mds/modrefine/stderr.txt
#PBS -N GalRefine
#PBS -M kevin.esoh@students.jkuat.ac.ke
#PBS -m b

cd /mnt/lustre/groups/CBBI1243/KEVIN/git/mds/modrefine/

GalaxyRefine -p $inpdb -t $pname -o 10

#bash galaxRefine.sh -n ocrlWt -p ../OCRL.B99990019.pdb
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
