#!/bin/bash

#SBATCH -J Homme		# job name
#SBATCH -o Homme.o%j		# output and error file name (%j expands to jobID)
#SBATCH -N @NODES@ -n @TASKS@	# total number of nodes and mpi tasks requested
#SBATCH -p normal-2mic		# queue (partition) -- normal, development, etc.
#SBATCH -t 00:30:00		# run time (hh:mm:ss) - 0.5 hours

ulimit -s unlimited

export OMP_NUM_THREADS=@HNT@
export MIC_OMP_NUM_THREADS=@MNT@
export MIC_PPN=@MPN@

[ -x /opt/intel/mic/bin/micinfo ] && /opt/intel/mic/bin/micinfo &> ./micinfo.log

