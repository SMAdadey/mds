#!/bin/bash

# Get input
gromMDS() {

   res=$1
   fle=$2
   f=$(basename $2)           # input file
   ff=$3		# force field
   nmds=$4		# Number of MD steps: 50000000 ; 2 * 50000000 = 100000 ps (100 ns)
   mdel=$5
   #fe="$(echo ${f##*.})"   # file extension
   #
   #if [[ $# != 1 ]]; then
   #   echo """
   #   Usage: gromMDS <pdb-file>
   #   """
   #elif [[ $# == 1 && $fe != "pdb" ]]; then
   #   echo "The input file is not a PDB file! Your file must end with .pdb"
   #
   #else
   
   # Initialize pbd input & output/intermediate filenames
   fb="${f/.pdb/}" # input base
   fc="${f/.pdb/_clean.pdb}" # clean file - water removed
   fp="${f/.pdb/_processed.gro}" # processed clean file
   fn="${f/.pdb/_newbox.gro}" # box file
   fs="${f/.pdb/_solv.gro}" # solvated file
   fsi="${f/.pdb/_solv_ions.gro}" # solvated file with ions added
   Epe="${f/.pdb/_potential.xvg}"; Epet="${Epe/.xvg/.txt}" # poten ener out files
   Et="${f/.pdb/_temperature.xvg}"; Ett="${Et/.xvg/.txt}" # temp out files
   Epr="${f/.pdb/_pressure.xvg}"; Eprt="${Epr/.xvg/.txt}" # press out files
   Ed="${f/.pdb/_density.xvg}"; Edt="${Ed/.xvg/.txt}" # density out files
   cpt="$fb.cpt" # checkpoint file
   wd="$(pwd)/${fb}_${ff}/"
   # Functions for generating parameter files
   writeIonMdp() {
   	echo """
   ; ions.mdp - used as input into grompp to generate ions.tpr
   ; Parameters describing what to do, when to stop and what to save
   integrator  = steep         ; Algorithm (steep = steepest descent minimization)
   emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
   emstep      = 0.01          ; Minimization step size
   nsteps      = 50000         ; Maximum number of (minimization) steps to perform
   
   ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
   nstlist         = 1         ; Frequency to update the neighbor list and long range forces
   cutoff-scheme	= Verlet    ; Buffered neighbor searching
   ns_type         = grid      ; Method to determine neighbor list (simple, grid)
   coulombtype     = cutoff    ; Treatment of long range electrostatic interactions
   rcoulomb        = 1.0       ; Short-range electrostatic cut-off
   rvdw            = 1.0       ; Short-range Van der Waals cut-off
   pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
   	"""
   }
   
   writeMinimMdp() {
   	echo """
   ; minim.mdp - used as input into grompp to generate em.tpr
   ; Parameters describing what to do, when to stop and what to save
   integrator  = steep         ; Algorithm (steep = steepest descent minimization)
   emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
   emstep      = 0.01          ; Minimization step size
   nsteps      = 50000         ; Maximum number of (minimization) steps to perform
   
   ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
   nstlist         = 1         ; Frequency to update the neighbor list and long range forces
   cutoff-scheme   = Verlet    ; Buffered neighbor searching
   ns_type         = grid      ; Method to determine neighbor list (simple, grid)
   coulombtype     = PME       ; Treatment of long range electrostatic interactions
   rcoulomb        = 1.0       ; Short-range electrostatic cut-off
   rvdw            = 1.0       ; Short-range Van der Waals cut-off
   pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
   	"""
   }
   
   writeNvtMdp() {
   echo -e """
   title                   = ${ff} $fb NVT equilibration
   define                  = -DPOSRES  ; position restrain the protein
   ; Run parameters
   integrator              = md        ; leap-frog integrator
   nsteps                  = 50000     ; 2 * 50000 = 100 ps
   dt                      = 0.002     ; 2 fs
   ; Output control
   nstxout                 = 500       ; save coordinates every 1.0 ps
   nstvout                 = 500       ; save velocities every 1.0 ps
   nstenergy               = 500       ; save energies every 1.0 ps
   nstlog                  = 500       ; update log file every 1.0 ps
   ; Bond parameters
   continuation            = no        ; first dynamics run
   constraint_algorithm    = lincs     ; holonomic constraints
   constraints             = h-bonds   ; bonds involving H are constrained
   lincs_iter              = 1         ; accuracy of LINCS
   lincs_order             = 4         ; also related to accuracy
   ; Nonbonded settings
   cutoff-scheme           = Verlet    ; Buffered neighbor searching
   ns_type                 = grid      ; search neighboring grid cells
   nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet
   rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
   rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
   DispCorr                = EnerPres  ; account for cut-off vdW scheme
   ; Electrostatics
   coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
   pme_order               = 4         ; cubic interpolation
   fourierspacing          = 0.16      ; grid spacing for FFT
   ; Temperature coupling is on
   tcoupl                  = V-rescale             ; modified Berendsen thermostat
   tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
   tau_t                   = 0.1     0.1           ; time constant, in ps
   ref_t                   = 300     300           ; reference temperature, one for each group, in K
   ; Pressure coupling is off
   pcoupl                  = no        ; no pressure coupling in NVT
   ; Periodic boundary conditions
   pbc                     = xyz       ; 3-D PBC
   ; Velocity generation
   gen_vel                 = yes       ; assign velocities from Maxwell distribution
   gen_temp                = 300       ; temperature for Maxwell distribution
   gen_seed                = -1        ; generate a random seed
   	"""
   }
   
   writeNptMdp() {
   	echo -e """
   title                   = ${ff} $fb NPT equilibration
   define                  = -DPOSRES  ; position restrain the protein
   ; Run parameters
   integrator              = md        ; leap-frog integrator
   nsteps                  = 50000     ; 2 * 50000 = 100 ps
   dt                      = 0.002     ; 2 fs
   ; Output control
   nstxout                 = 500       ; save coordinates every 1.0 ps
   nstvout                 = 500       ; save velocities every 1.0 ps
   nstenergy               = 500       ; save energies every 1.0 ps
   nstlog                  = 500       ; update log file every 1.0 ps
   ; Bond parameters
   continuation            = yes       ; Restarting after NVT
   constraint_algorithm    = lincs     ; holonomic constraints
   constraints             = h-bonds   ; bonds involving H are constrained
   lincs_iter              = 1         ; accuracy of LINCS
   lincs_order             = 4         ; also related to accuracy
   ; Nonbonded settings
   cutoff-scheme           = Verlet    ; Buffered neighbor searching
   ns_type                 = grid      ; search neighboring grid cells
   nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
   rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
   rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
   DispCorr                = EnerPres  ; account for cut-off vdW scheme
   ; Electrostatics
   coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
   pme_order               = 4         ; cubic interpolation
   fourierspacing          = 0.16      ; grid spacing for FFT
   ; Temperature coupling is on
   tcoupl                  = V-rescale             ; modified Berendsen thermostat
   tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
   tau_t                   = 0.1     0.1           ; time constant, in ps
   ref_t                   = 300     300           ; reference temperature, one for each group, in K
   ; Pressure coupling is on
   pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
   pcoupltype              = isotropic             ; uniform scaling of box vectors
   tau_p                   = 2.0                   ; time constant, in ps
   ref_p                   = 1.0                   ; reference pressure, in bar
   compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
   refcoord_scaling        = com
   ; Periodic boundary conditions
   pbc                     = xyz       ; 3-D PBC
   ; Velocity generation
   gen_vel                 = no        ; Velocity generation is off
   	"""
   }
   
   writeMdMdp() {
   	echo -e """
   title                   = ${ff} $fb NPT equilibration
   ; Run parameters
   integrator              = md        ; leap-frog integrator
   nsteps                  = ${nmds}     ; 2 * 50000000 = 100000 ps (100 ns)
   dt                      = 0.002     ; 2 fs
   ; Output control
   nstxout                 = 0         ; suppress bulky .trr file by specifying
   nstvout                 = 0         ; 0 for output frequency of nstxout,
   nstfout                 = 0         ; nstvout, and nstfout
   nstenergy               = 5000      ; save energies every 10.0 ps
   nstlog                  = 5000      ; update log file every 10.0 ps
   nstxout-compressed      = 5000      ; save compressed coordinates every 10.0 ps
   compressed-x-grps       = System    ; save the whole system
   ; Bond parameters
   continuation            = yes       ; Restarting after NPT
   constraint_algorithm    = lincs     ; holonomic constraints
   constraints             = h-bonds   ; bonds involving H are constrained
   lincs_iter              = 1         ; accuracy of LINCS
   lincs_order             = 4         ; also related to accuracy
   ; Neighborsearching
   cutoff-scheme           = Verlet    ; Buffered neighbor searching
   ns_type                 = grid      ; search neighboring grid cells
   nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
   rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
   rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
   ; Electrostatics
   coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
   pme_order               = 4         ; cubic interpolation
   fourierspacing          = 0.16      ; grid spacing for FFT
   ; Temperature coupling is on
   tcoupl                  = V-rescale             ; modified Berendsen thermostat
   tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
   tau_t                   = 0.1     0.1           ; time constant, in ps
   ref_t                   = 300     300           ; reference temperature, one for each group, in K
   ; Pressure coupling is on
   pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
   pcoupltype              = isotropic             ; uniform scaling of box vectors
   tau_p                   = 2.0                   ; time constant, in ps
   ref_p                   = 1.0                   ; reference pressure, in bar
   compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
   ; Periodic boundary conditions
   pbc                     = xyz       ; 3-D PBC
   ; Dispersion correction
   DispCorr                = EnerPres  ; account for cut-off vdW scheme
   ; Velocity generation
   gen_vel                 = no        ; Velocity generation is off
   	"""
   }
  
   fe="$(echo ${f##*.})"   # file extension

function make_params() {
   # Generate Molecular Dynamic Parameter Files
   writeIonMdp > ions.mdp # for solvation
   writeMinimMdp > minim.mdp # for energy minimization
   writeNvtMdp > nvt.mdp # for NVT equilibration
   writeNptMdp > npt.mdp # for NPT equilibration
   writeMdMdp > md.mdp # for production MD
}

function qsub_prep() {
echo -e """#!/bin/bash
#PBS -l select=1:ncpus=24:mpiprocs=24
#PBS -l walltime=96:00:00
#PBS -q smp
#PBS -P CBBI1243
#PBS -o ${wd}stdout.txt
#PBS -e ${wd}stderr.txt
#PBS -m b
#PBS -N Gromacs_$fb
#PBS -M kevin.esoh@students.jkuat.ac.ke

#MODULEPATH=/opt/gridware/bioinformatics/modules:\$MODULEPATH
#source /etc/profile.d/modules.sh

### # module add/load # ###
module add chpc/BIOMODULES
module load gromacs/5.1.4-openmpi_1.10.2-intel16.0.1
module load R/3.6.0-gcc7.2.0

#OMP_NUM_THREADS=1 #turn off OpenMP (also -ntomp on commandline)

#-nt -ntmpi number of threads and number of mpi threads

NP=\`cat \${PBS_NODEFILE} | wc -l\`
#NP=24
echo \"Number of Processes: \$NP\"

mdr=\"gmx_mpi mdrun\"
#ARGS=\"-s X -deffnm Y\"

cd ${wd}
#mpirun -np \${NP} -machinefile \${PBS_NODEFILE} \${mdr} \${ARGS}""" > qsub_prep
}

function nem_md() {
echo -e """#!/usr/bin/env bash

# Remove water molecules (crystal)
grep -v \"HOH\" $f > $fc

# Create required files; topology, position restraint, post-processed structure
gmx pdb2gmx -f $fc -o $fp -water $mdel -ignh -ff $ff	 # AMBER99SB ff #ignh: ignore H atoms in PDB file

# Define unit cell & add solvent
gmx editconf -f $fp -o $fn -c -d 1.0 -bt dodecahedron
gmx solvate -cp $fn -cs spc216.gro -o $fs -p topol.top
gmx grompp -f ions.mdp -c $fs -p topol.top -o ions.tpr
echo 13 | gmx genion -s ions.tpr -o $fsi -p topol.top -pname NA -nname CL -neutral

# Energy minimization
gmx grompp -f minim.mdp -c $fsi -p topol.top -o em.tpr
mpirun -np 1 gmx mdrun -v -deffnm em
echo -e \"Potential\\\n0\"  | gmx energy -f em.edr -o $Epe
grep -v -e \"#\" -e \"@\" $Epe > $Epet

# Equilibration Phase 1: NVT (Energy and Temperature)
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
mpirun -np 1 gmx mdrun -v -deffnm nvt
echo -e \"Temperature\\\n0\" | gmx energy -f nvt.edr -o $Et
grep -v -e \"#\" -e \"@\" $Et > $Ett

# Equilibration Phase 2: NPT (Pressure and Density)
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
mpirun -np 1 gmx mdrun -v -deffnm npt
echo -e \"Pressure\\\n0\" | gmx energy -f npt.edr -o $Epr
grep -v -e \"#\" -e \"@\" $Epr > $Eprt
echo -e \"Density\\\n0\" | gmx energy -f npt.edr -o $Ed
grep -v -e \"#\" -e \"@\" $Ed > $Edt

# Run Production MD
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
mpirun -np 1 gmx mdrun -v -deffnm md_0_1""" > nem_md
}

function em_md() {
echo -e """
# Remove water molecules (crystal)
grep -v \"HOH\" $f > $fc

# Create required files; topology, position restraint, post-processed structure
gmx_mpi pdb2gmx -f $fc -o $fp -water $mdel -ignh -ff $ff  # AMBER99SB ff

# Define unit cell & add solvent
gmx_mpi editconf -f $fp -o $fn -c -d 1.0 -bt dodecahedron
gmx_mpi solvate -cp $fn -cs spc216.gro -o $fs -p topol.top
gmx_mpi grompp -f ions.mdp -c $fs -p topol.top -o ions.tpr
echo 13 | gmx_mpi genion -s ions.tpr -o $fsi -p topol.top -pname NA -nname CL -neutral

# Energy minimization
gmx_mpi grompp -f minim.mdp -c $fsi -p topol.top -o em.tpr
time \${mdr} -v -cpi -maxh 95 -ntomp \${NP} -deffnm em #gmx mdrun
echo -e \"Potential\\\n0\" | gmx_mpi energy -f em.edr -o $Epe
grep -v -e \"#\" -e \"@\" $Epe > $Epet

# Equilibration Phase 1: NVT (Energy and Temperature)
gmx_mpi grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
time \${mdr} -v -maxh 95 -ntomp \${NP} -deffnm nvt #gmx mdrun
echo -e \"Temperature\\\n0\" | gmx_mpi energy -f nvt.edr -o $Et
grep -v -e \"#\" -e \"@\" $Et > $Ett

# Equilibration Phase 2: NPT (Pressure and Density)
gmx_mpi grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
time \${mdr} -v -maxh 95 -ntomp \${NP} -deffnm npt #gmx mdrun
echo -e \"Pressure\\\n0\" | gmx_mpi energy -f npt.edr -o $Epr
grep -v -e \"#\" -e \"@\" $Epr > $Eprt
echo -e \"Density\\\n0\" | gmx_mpi energy -f npt.edr -o $Ed
grep -v -e \"#\" -e \"@\" $Ed > $Edt

# Run Production MD
gmx_mpi grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
time \${mdr} -v -maxh 95 -ntomp \${NP} -deffnm md_0_1 #gmx mdrun """ > em_md
}

function analysis() {
echo -e """
# Analysis
echo 1 0 | gmx_mpi trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center
echo 4 4 | gmx_mpi rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns
grep -v -e \"#\" -e \"@\" rmsd.xvg > rmsd.txt

echo 4 4 | gmx_mpi rms -s em.tpr -f md_0_1_noPBC.xtc -o rmsd_xtal.xvg -tu ns
grep -v -e \"#\" -e \"@\" rmsd_xtal.xvg > rmsd_xtal.txt

echo 1 | gmx_mpi gyrate -s md_0_1.tpr -f md_0_1_noPBC.xtc -o gyrate.xvg
grep -v -e \"#\" -e \"@\" gyrate.xvg > gyrate.txt

#gmx_mpi rama -f em.gro -s em.tpr -o ramachan.xvg # Ramachandran Plot for crystal struct""" > analysis
}

function plot() {
echo -e """
# Generate plots
Rscript ${HOME}/Git/mds/mds_plot.R $Epet $Ett $Eprt $Edt rmsd.txt rmsd_xtal.txt gyrate.txt gyrate.txt""" > plot
}

function nhrestart() {
echo -e """
# Restart Production MD
mpirun -np 1 gmx mdrun -v -cpi md_0_1""" > nhrestart
}

function hrestart() {
echo -e """
# Restart Production MD
time \${mdr} -s md_0_1 -v -maxh 95 -ntomp \${NP} -cpi md_0_1 #gmx mdrun""" > hrestart
}

function msg() {
      sleep 1
      if [[ ($res == "hpc") || ($res == "hrestart") ]]; then
          echo -e "\nMolecular Dynamics Parameter (mdp) files created in your current directory!"
          sleep 1
          echo -e "Gromacs job created: '${fb}.${ff}.${res}.qsub' Please submit with 'qsub ${fb}.${ff}.${res}.qsub'\n"
      else
          echo -e "\nMolecular Dynamics Parameter (mdp) files created in ${fb}_${ff}"
          sleep 1
          echo -e "Gromacs job created: '${fb}.${ff}.${res}.sh' currently in ${fb}_${ff}/ Please run with './${fb}.${ff}.${res}.sh'\n"
      fi
}

function prep_all() { 
      qsub_prep; nem_md; em_md; nhrestart; hrestart; analysis; plot
}

function rm_all() {
      rm qsub_prep em_md nem_md nhrestart hrestart analysis plot
}

   if [[ $# != 5 ]]; then
      echo """
      Usage: gromMDS <[hpc | hrestart | nhpc | nrestart]> <pdb-file> <force-field> <nsteps> <water-model>
             
             hpc                : Run MDS on HPC
             hrestart           : Restart MDS on HPC
             nhpc               : Run MDS on local terminal or interactively
             nhrestart          : Restart MDS interactively
	     
	     pdb-file		: The coordinate file to run MDS for
	     force-field	: One of several force fields ported to GROMACS

	     ----------------------------------------------------------------------------------------------------------
		 1: AMBER03 protein, nucleic AMBER94 (Duan et al., J. Comp. Chem. 24, 1999-2012, 2003)
		 2: AMBER94 force field (Cornell et al., JACS 117, 5179-5197, 1995)
		 3: AMBER96 protein, nucleic AMBER94 (Kollman et al., Acc. Chem. Res. 29, 461-469, 1996)
		 4: AMBER99 protein, nucleic AMBER94 (Wang et al., J. Comp. Chem. 21, 1049-1074, 2000)
		 5: AMBER99SB protein, nucleic AMBER94 (Hornak et al., Proteins 65, 712-725, 2006)
		 6: AMBER99SB-ILDN protein, nucleic AMBER94 (Lindorff-Larsen et al., Proteins 78, 1950-58, 2010)
		 7: AMBERGS force field (Garcia & Sanbonmatsu, PNAS 99, 2782-2787, 2002)
		 8: CHARMM27 all-atom force field (CHARM22 plus CMAP for proteins)
		 9: CHARMM36 all-atom force field (March 2019)
		10: GROMOS96 43a1 force field
		11: GROMOS96 43a2 force field (improved alkane dihedrals)
		12: GROMOS96 45a3 force field (Schuler JCC 2001 22 1205)
		13: GROMOS96 53a5 force field (JCC 2004 vol 25 pag 1656)
		14: GROMOS96 53a6 force field (JCC 2004 vol 25 pag 1656)
		15: GROMOS96 54a7 force field (Eur. Biophys. J. (2011), 40,, 843-856, DOI: 10.1007/s00249-011-0700-9)
		16: OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)
		
		NB: Case sensitive. Must enter lower case on commandline not upper case.
		Example:
		  - OPLS-AA/L = oplsaa
		  - AMBER99SB = amber99sb
             ----------------------------------------------------------------------------------------------------------

	     nsteps		: Number of MD steps to run
	     			  50000000 (50 million) = 2 * 50000000 = 100000 ps = 100 ns
	     water-model	: Water model to use: select, none, spc, spce, tip3p, tip4p, tip5p, tips3p

      """
   elif [[ $# == 5 && $fe != "pdb" ]]; then
      echo -e "\nThe input file is not a PDB file! Your file must end with .pdb\n"

   elif [[ $# == 5 && $res == "nhpc" && $fe == "pdb" ]]; then
	mkdir -p ${fb}_${ff}; cp $fle ${wd}/; cd $wd
        make_params; prep_all
        cat nem_md analysis plot | sed 's/gmx_mpi/gmx/g' > ${fb}.${ff}.${res}.sh
        chmod 755 ${fb}.${ff}.${res}.sh
	msg; rm_all
	#./${fb}.${ff}.${res}.sh
   elif [[ $# == 5 && $res == "nhrestart" && $fe == "pdb" ]]; then
        mkdir -p ${fb}_${ff}; cp $fle ${wd}/; cd $wd
        make_params; prep_all
        echo -e "#!/usr/bin/env bash\nmdr=\"mpirun -np 2 gmx mdrun\"" > ${fb}.${ff}.${res}.sh
        cat nhrestart analysis plot | sed 's/gmx_mpi/gmx/g' >> ${fb}.${ff}.${res}.sh
        chmod 755 ${fb}.${ff}.${res}.sh
	msg; rm_all
        #./${fb}.${ff}.${res}.sh
   elif [[ $# == 5 && $res == "hpc" && $fe == "pdb" ]]; then
	mkdir -p ${fb}_${ff}; cp $fle ${wd}/; cd ${wd}
        make_params; prep_all
        cat qsub_prep em_md analysis plot > ${fb}.${ff}.${res}.qsub
	msg; rm_all
        #qsub ${fb}.${ff}.${res}.qsub
   elif [[ $# == 5 && $res == "hrestart" && $fe == "pdb" ]]; then
        mkdir -p ${fb}_${ff}; cp $fle ${wd}/; cd ${wd}
        make_params; prep_all
        cat qsub_prep hrestart analysis plot > ${fb}.${ff}.${res}.qsub
	msg; rm_all
        #qsub ${fb}.${ff}.${res}.qsub
   fi

}
