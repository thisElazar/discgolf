#!/bin/sh
#SBATCH --job-name=anchors
#SBATCH --partition=general
#SBATCH --array=0-20%20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --output=anchor_%a.log
# 3 anchor discs x 7 angles = 21 solves on the validated coarse mesh.
~/cfd/bin/micromamba run -p ~/cfd/env sh "$HOME/cfd/disc/cfd/anchors/run_anchor.sh"
