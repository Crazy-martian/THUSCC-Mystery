#!/bin/bash

#SBATCH -J MYSTERY-RUN
#SBATCH -p cnall        
#SBATCH -N 1          
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=56
#SBATCH -o logs/stdout.%j    
#SBATCH -e logs/stderr.%j 

source ~/.bashrc
conda activate mystery_env


echo "X Ray Image task: "


echo "================================================================"
echo "[$(date +"%Y-%m-%d %H:%M:%S.%3N")] Running..."
echo "================================================================"

# python ratio_calculate.py

export LD_LIBRARY_PATH=$CFITSIO_INSTALL_DIR/lib:$LD_LIBRARY_PATH

./ratio_calc

echo "================================================================"
echo "[$(date +"%Y-%m-%d %H:%M:%S.%3N")] Checking..."
echo "================================================================"

python score.py

echo "================================================================"
echo "[$(date +"%Y-%m-%d %H:%M:%S.%3N")] Done."
echo "================================================================"
