#!/bin/bash
# Mesh the three anchor discs on the HEAD node (same settings as the bigsweep:
# shared template system/ = validated coarse mesh, y+~30 wall-function regime).
# Run:  ~/cfd/bin/micromamba run -p ~/cfd/env bash mesh_anchors.sh
BASE="$HOME/cfd/disc/cfd"; AN="$BASE/anchors"; TPL="$BASE/bigsweep/template"
export AN TPL

ls -d "$AN"/geoms/anchor_* | xargs -P 3 -I{} bash -c '
  g="$1"; cd "$g" || exit 1
  ln -sfn "$TPL/system" system
  mkdir -p constant
  cp "$TPL/constant/transportProperties" "$TPL/constant/turbulenceProperties" constant/ 2>/dev/null
  blockMesh            > log.blockMesh 2>&1 || { echo "FAIL block  $(basename $g)"; exit 1; }
  snappyHexMesh -overwrite > log.snappy 2>&1 || { echo "FAIL snappy $(basename $g)"; exit 1; }
  echo "OK $(basename $g)"
' _ {}

for g in "$AN"/geoms/anchor_*; do
  [ -f "$g/constant/polyMesh/owner" ] && echo "meshed $(basename $g)" || echo "UNMESHED $(basename $g)"
done
