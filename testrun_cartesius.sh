#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH -J testrun

cd $HOME/GalapagosNEMO/scripts/particleruns
echo 'Executing program ...'

python galapagosrun_forward.py

echo 'Finished computation.'
