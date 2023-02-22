#!/bin/bash
#SBATCH --job-name=run
#SBATCH --output=output_log/run.txt
#SBATCH --mem=8G
#SBATCH --time=1000:00:00
#SBATCH -c 4
#SBATCH --nodelist=XXX

echo "Beginning of script"
date

python3 riskratio_tables.py

echo "End of script"
date