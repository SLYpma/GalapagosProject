#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH -J nemo2008

echo 'Executing program ...'

python galapagosrun_bwd_nemo_2008.py

echo 'Finished computation.'
