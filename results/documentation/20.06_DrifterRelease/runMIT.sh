#!/bin/bash
#SBATCH -t 2-00:00:00
#SBATCH -N 1
#SBATCH -J ReleaseDrifter_MITgcm
#SBATCH --mail-type=begin       
#SBATCH --mail-type=end         
#SBATCH --mail-type=fail        
#SBATCH --mail-user=s.l.ypma@uu.nl

echo 'Executing program ...'

python drifterrun_fwd_MITgcm4km.py

echo 'Finished computation.'
