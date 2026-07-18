#!/bin/sh
#SBATCH --job-name=discsweep
#SBATCH --partition=general
#SBATCH --array=0-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --output=%x_%a.log
# 21 independent angles fan out across mgd[1-5] (8 tasks/node).
~/cfd/bin/micromamba run -p ~/cfd/env sh "$HOME/cfd/disc/cfd/run_alpha.sh"
