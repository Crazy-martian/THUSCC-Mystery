#!/bin/bash

#SBATCH -J MYSTERY-OPTIMIZED
#SBATCH -p cnall         # 使用队列
#SBATCH -N 1             # 使用节点数
#SBATCH --ntasks-per-node=1  # 每个节点使用的核数
#SBATCH --cpus-per-task=56
#SBATCH -o logs/stdout.%j      # 标准输出
#SBATCH -e logs/stderr.%j      # 错误输出

source ~/.bashrc
conda activate mystery_env



echo "X Ray Image task:"


echo "================================================================"
echo "[$(date +"%Y-%m-%d %H:%M:%S.%3N")] Running..."
echo "================================================================"

python ratio_calculate.py

echo "================================================================"
echo "[$(date +"%Y-%m-%d %H:%M:%S.%3N")] Checking..."
echo "================================================================"

python score.py

echo "================================================================"
echo "[$(date +"%Y-%m-%d %H:%M:%S.%3N")] Done."
echo "================================================================"
