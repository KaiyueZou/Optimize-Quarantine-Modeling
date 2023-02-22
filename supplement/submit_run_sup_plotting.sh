#!/bin/bash
#SBATCH --job-name=run_sup_plotting
#SBATCH --output=output_run_sup_plotting.txt
#SBATCH --time=10000:00:00
#SBATCH --mem=8G
#SBATCH --nodelist=XXX
#SBATCH -c 4

echo "Beginning of script"
date
TAXDIR=/.../

cd $TAXDIR
python3 $TAXDIR/riskratio_plotting_clean_new.py

echo "End of script"
date