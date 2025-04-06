#!/bin/zsh
#SBATCH --job-name=MD_simulation
#SBATCH --partition=mixed
#SBATCH --nodelist=DL380-2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=30-00:00:00
#SBATCH --output=file_1.out
#SBATCH --error=file_1.err
#SBATCH --mem=60G


## 环境配置方法
## 1 直接导入环境变量

## 2 使用environment-module加载  %% 需要提前准备好环境
. ~/load.sh
source activate yrydmd
## 3 在提交作业之前，用shell激活环境 （以下代码在shell中执行，不用写入sbatch脚本

### 作业命令
## 1. 执行并行计算

python3 /home/lihuaxiang/pysimm/PyPA/main.py

# mpirun -np 32 lmp_mpi -sf gpu -pk gpu 1 -in run.in.npt
