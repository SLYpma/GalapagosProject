#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -N 1
#SBATCH -J 4km_wstokes

echo 'Executing program ...'

python galapagosrun_bwd_4km_2008.py

echo 'Finished computation.'
