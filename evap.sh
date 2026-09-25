#!/bin/bash
#SBATCH -p d8168         # Queue
#SBATCH -N 1          # Node count required for the job
#SBATCH -n 48           # Number of tasks to be launched

# USE:      ./evap.sh <cutoff> <nmax> [start]
# USE:      sbatch evap.sh <cutoff> <nmax> [start]  (SLURM)
#  
# <cutoff> = cutoff distance for evaporation
# <nmax> = maximal step number for evaporation
# [start]     = the *optional* third argument "start" tells the script whether you are
#                 starting an evaporation; if no third argument is given to the script,
#                 the evaporation will be a restart of a unfinished evaporation 

# List of files that need to be in this folder
# COORDINATES: out.gro
# TOPOLOGY: topol.top, P.itp, S.itp
# RUN PARAMETERS: em.mdp, md.mdp
# OTHERS: out.tpr; sum.txt

# Load GROMACS
source /public1/AMD/app/module/5.0.1/init/sh-intel
module load gromacs/4.5.5-intel2016

# System-dependent parameters
m0=$(expr 6 \* $(grep "mR      Rm" out.gro | wc -l)) # number of CG beads for the 3D ladder CANAL polymers in the system

# Check inputs
cutoff="$1"		# It has an unit of nm
nmax="$2"		# it should be a positive integer
size_1=${#cutoff}		
size_2=${#nmax}		

if [ $size_1 -eq 0 -o $size_2 -eq 0 ] ; then
  echo ""
  echo "Missing cutoff distance or maximal step number. Check that"
  exit
fi

if [ $nmax -gt 0 -a $nmax -le 100 ]; then
  nmax=$(printf "%03d" $nmax)   # makes sure also "$nmax" is a 3 digits long variable
else
  echo "The maximal step number is too big for the computational effort or too small to be reasonable. Change it"
  exit
fi

# START a new evaporation OR RESTART an unfinished evaporation?
is_a_start="$3"
size_3=${#is_a_start}

if [ "$is_a_start" == "start" ] ; then
####################
## It's a start! ##
####################
  rm -rf step*
  mkdir step000
  cd step000
  cp ../out.gro out.gro
  cp ../out.tpr out.tpr
  cp ../index.ndx index.ndx
  cd -
# Define where to start the evaporation loop from (i.e., step000, since it's a start)
  starting_i=000

elif [ $size_3 = 0 ] ; then
######################
## It's a restart! ##
######################
  echo -e "\n That's an evaporation restart"
# Where do we have to start from? Look for folders and the file 'out.gro'
  for i in {001..100} ; do

    find ./step$i -name 'out.gro' 1> /dev/null
    findres=$?

    if [ $findres != 0 ] ; then
# No folder has been found for this step. The previous step had completed.
      echo "Folder not found; we have to restart from" $i
# Define where to start the evaporation loop from
      previous=$(expr $i - 1)               # previous step variable
      previous=$(printf "%03d" $previous)   # makes sure also "previous" is a 3 digits long variable 
      echo "Let's move to the previous folder (number" $previous ")"
      starting_i=$previous
      break
    fi

    found_a_outfile=$(find ./step$i -name 'out.gro' | wc -l)

    if [ $found_a_outfile = 0 ] ; then
   # No file 'finished' has been found in the current folder --> start from CHECKPOINT (state_prev.cpt)
      echo "we have to restart from step" $i
      # Define where to start the evaporation loop from
      echo "Let's move to the folder (number" $i ")"
      # Let's treat "i" as the previous step, for compatibility with EVAPORATION
      previous=$(expr $i - 1)               # previous step variable
      previous=$(printf "%03d" $previous)   # makes sure also "previous" is a 3 digits long variable 
      starting_i=$previous
      break
    elif [ $found_a_outfile = 1 ] ; then
      echo "step" $i "had finished"
    else
      echo "Something is wrong? More than one 'finished' file has been found in folder 'step"$i"'"
      exit
    fi

  done

else
######################
## Error! ##
######################
   echo "If you want to start a new evaporation process, pass "start" (as third argument) to the script."
   echo "If you want to restart a unfinished evaporation process, no third argument is required."
   exit
fi

######################
## EVAPORATION ##
######################
i=$starting_i

while [ "$i" -lt "$nmax" ] ; do

  found_a_outfile=$(find ./step$i -name 'out.gro' | wc -l)

  if [ $found_a_outfile = 1 ] ; then
    echo "step" $i "had finished"
  else
    echo "There's something wrong. The evaporation stops at 'step"$i"'"
    break
  fi

  j=$(expr $i + 1)        # "j" = 'next step'
  j=$(printf "%03d" $j)   # makes sure also "j" is a 3 digits long variable

  if [[ -d "step$j" ]]; then
    echo step$j "already exists"
  else
    mkdir step$j
  fi

  cd step$j

# obtain number of CG beads (including P and S) in the system to evaporate
  n1=$(sed -n '2p' "../step$i/out.gro" | tr -d ' ')

# generate the index file for the CG beads in the remaining system after evaporation
  g_select_mpi -f ../step$i/out.gro -s ../step$i/out.tpr -n ../step$i/index.ndx -select 'group "P" or (resname S and within '$cutoff' of group "P")' -on tmp.ndx

# obtain the coordinate (.gro) file for the remaining system after evaporation
  echo -e "0\n" | \trjconv_mpi -f ../step$i/out.gro -s ../step$i/out.tpr -n tmp.ndx -o in.gro

# judge if evaporation does occur
  n2=$(sed -n '2p' "in.gro" | tr -d ' ')
  if [ $(expr $n1 - $n2) -eq 0 ]; then
    echo ""
    echo "No more solvents can be evaporated. You can need to perform additional relaxations or decrease the cutoff distance."
    exit
  fi

# complete the forcefield (.top) file for the remaining system after evaporation
  cp ../topol.top topol.top
  cp ../P.itp P.itp
  cp ../S.itp S.itp

# make the index (.ndx) file for the goups named P and S that are mentioned in the parameter (.mdp) files 
   m1=$(expr $n2 - $m0)

  if [ $m1 != 0 ]; then
    echo -e "2 | 3 | 4 | 5 | 6 | 7\nname 9 P\nq" | \make_ndx_mpi -f in.gro -o
    echo "S             $m1" >> topol.top
  else
    echo -e "2 | 3 | 4 | 5 | 6 | 7\nname 8 P\nq" | \make_ndx_mpi -f in.gro -o
  fi

# prepare and execute steepest descent (sd) and conjugate gradient (cg) energy minimizations (em)
  grompp_mpi -f ../em.mdp -p -c in.gro -n -o in_em.tpr -po in_em.mdp -maxwarn 10
  wait
  mpirun -np 1 mdrun_mpi -s in_em.tpr -deffnm in_em -v

# prepare and execute molecular dynamics (md) relaxation simulation in constant-NVT ensemble
  grompp_mpi -f ../md.mdp -p -c in_em.gro -n -o out.tpr -po out.mdp -maxwarn 10
  wait
  mpirun -np 4 mdrun_mpi -s out.tpr -deffnm out -v -pd

  echo "Run step number" $j
  cd -

  echo "$j	$m1" >> sum.txt

  i=$(expr $i + 1)
  i=$(printf "%03d" $i)   # makes sure the new "i" is a 3 digits long variable 
done

echo "Evaporation has finished. Please check if it is normal"