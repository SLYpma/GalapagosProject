#!/bin/bash
# SGE Options
#$ -S /bin/bash
# Shell environment forwarding
#$ -V
# Job Name
#$ -N galapagos_forward_test
# Notifications
#$ -M s.l.ypma@uu.nl
# When notified
#$ -m es
# Set memory limit
#$ -l h_vmem=20G
# Set runtime limit
#$ -l h_rt=20:00:00
# run the job on the queue for long-running processes:
#$ -q all.q

echo 'Executing program ...'

python /scratch/SLYpma/GalapagosNEMO/scripts/particleruns/galapagosrun_forward.py

echo 'Finished computation.'
