#!/bin/sh
#SBATCH --job-name=bigsweep
#SBATCH --partition=general
#SBATCH --array=0-671%20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --output=bigsweep_%a.log
# 96 geometries x 7 angles = 672 solves. %20 caps to 20 concurrent (matches
# the ~4 physical cores x 5 nodes), so it never oversubscribes and coexists
# with anything else queued. Submit now; it waits behind the validation job.
~/cfd/bin/micromamba run -p ~/cfd/env sh "$HOME/cfd/disc/cfd/bigsweep/run_case.sh"
