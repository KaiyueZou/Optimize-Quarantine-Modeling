#!/bin/bash
#SBATCH --job-name=run_sup
#SBATCH --output=output_run_sup.txt
#SBATCH --time=10000:00:00
#SBATCH --mem=8G
#SBATCH --nodelist=XXX
#SBATCH -c 4

echo "Beginning of script"
date
TAXDIR=/.../

cd $TAXDIR
python3 $TAXDIR/riskratio_tables.py

echo "End of script"
date