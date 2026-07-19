#!/bin/bash
# Mesh every big-sweep geometry on the HEAD NODE (110 GB, 24 cores) — the
# compute nodes have too little RAM for 0.5M-cell snappy. Parallelised with
# xargs; each geometry meshes serially into its own constant/polyMesh.
# Run:  ~/cfd/bin/micromamba run -p ~/cfd/env bash mesh_all.sh [parallelism]
BASE="$HOME/cfd/disc/cfd"; BS="$BASE/bigsweep"; TPL="$BS/template"
J="${1:-16}"
export BASE BS TPL

ls -d "$BS"/geoms/geom_* | xargs -P "$J" -I{} bash -c '
  g="$1"; cd "$g" || exit 1
  [ -f constant/polyMesh/owner ] && { echo "SKIP (meshed) $(basename $g)"; exit 0; }
  ln -sfn "$TPL/system" system
  mkdir -p constant
  cp "$TPL/constant/transportProperties" "$TPL/constant/turbulenceProperties" constant/ 2>/dev/null
  blockMesh            > log.blockMesh 2>&1 || { echo "FAIL block  $(basename $g)"; exit 1; }
  snappyHexMesh -overwrite > log.snappy 2>&1 || { echo "FAIL snappy $(basename $g)"; exit 1; }
  echo "OK $(basename $g) $(grep -m1 "cells:" log.snappy 2>/dev/null | tr -s " ")"
' _ {}

echo "=== summary ==="
ok=0; miss=0
for g in "$BS"/geoms/geom_*; do
  if [ -f "$g/constant/polyMesh/owner" ]; then ok=$((ok+1)); else miss=$((miss+1)); echo "  UNMESHED $(basename $g)"; fi
done
echo "meshed: $ok   unmeshed: $miss"
