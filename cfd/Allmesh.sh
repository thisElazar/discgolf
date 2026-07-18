#!/bin/sh
# Build the mesh once (blockMesh background + snappyHexMesh around disc.stl).
# Run on a compute node inside the OpenFOAM env, e.g.:
#   srun -N1 -n1 -c8 -p general ~/cfd/bin/micromamba run -p ~/cfd/env sh Allmesh.sh
cd "$(dirname "$0")/case" || exit 1
rm -rf processor* [1-9]* 0.[0-9]* log.* constant/polyMesh

blockMesh            > log.blockMesh   2>&1 || { echo "blockMesh FAILED";   tail -15 log.blockMesh;   exit 1; }
snappyHexMesh -overwrite > log.snappy 2>&1 || { echo "snappyHexMesh FAILED"; tail -25 log.snappy;      exit 1; }
checkMesh -constant  > log.checkMesh  2>&1

echo "=== snappy tail ==="; tail -4 log.snappy
echo "=== checkMesh ==="
grep -E "cells:|faces:|points:|Max aspect|Max non-ortho|Max skewness|Mesh OK|\*\*\*" log.checkMesh
