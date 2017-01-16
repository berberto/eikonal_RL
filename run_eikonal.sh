#!/bin/bash
#TO BE specified only in case of usage of a given queue (now commented)
##PBS -q reserved2
##PBS -q long
#PBS -q regular

#define the job array
##PBS -t 0-19

#defne the number of nodes (and cores)
#PBS -l nodes=1:ppn=20

#define the walltime
##PBS -l walltime=24:00:00
#PBS -l walltime=12:00:00

#define the error output file as JOBNAME.e[JOBID] || JOBNAME.o[JOBID]
#PBS -j oe

#send me an e-mail when job ends
#PBS -m abe
#PBS -M apezzotta@sissa.it

n='1000'

alpha='.1'
beta='.5'
eta='.1'

#Define the JOB name
#PBS -N $n_000-019_$alpha_$beta_$eta

#IT HAS TO BE HERE
#PBS -T flush_cache
#PBS
#

# set up environment to import pre-compiled executables from AK.
#module load intel

#set the current directory as the working directory (where the job is submitted from)
#by default the job will start on the $HOME directory on the remote node
PBS_O_WORK=/scratch/apezzotta/eikonal/main

cd $PBS_O_WORK

module load openmpi/1.8.3/intel/14.0

gvals=()
for i in `seq 0 9`; do
    gvals+=("0.00"$i"0")
done
for i in `seq 10 19`; do
    gvals+=("0.0"$i"0")
done

mpirun -np ${#gvals[@]} eikonal2d_mpi $n $alpha $beta $eta ${gvals[@]}
