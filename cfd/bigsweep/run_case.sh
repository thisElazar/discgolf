#!/bin/sh
# One Slurm array task = one (geometry, angle) pair.
#   task_id = geom_index * N_ALPHA + alpha_index
# Builds a work dir that SYMLINKS that geometry's pre-built mesh, stamps the
# freestream/force directions with THAT geometry's normalization (area, dia,
# CofR all differ per geometry), and runs simpleFoam serial.
BASE="$HOME/cfd/disc/cfd"; BS="$BASE/bigsweep"; TPL="$BS/template"
ALPHAS="-4 -1 2 5 8 11 14"; NA=7

# TASK_OFFSET lets chunked submissions cover task ids beyond Slurm's
# MaxArraySize (1001 here): real id = array index + offset.
T=$(( ${SLURM_ARRAY_TASK_ID:-0} + ${TASK_OFFSET:-0} ))
GI=$((T / NA)); AI=$((T % NA))
A=$(echo $ALPHAS | cut -d' ' -f $((AI + 1)))
GID=$(printf "geom_%03d" "$GI")
GEOM="$BS/geoms/$GID"

[ -f "$GEOM/constant/polyMesh/owner" ] || { echo "no mesh for $GID (skip)"; exit 0; }

WORK="$BS/solve/${GID}_a${A}"
rm -rf "$WORK"; mkdir -p "$WORK/constant"
cp -r "$TPL/0" "$WORK/0"
cp -r "$TPL/system" "$WORK/system"
cp "$TPL/constant/transportProperties" "$TPL/constant/turbulenceProperties" "$WORK/constant/"
ln -s "$GEOM/constant/polyMesh" "$WORK/constant/polyMesh"

cd "$BASE" || exit 1
python stamp_alpha.py --case "$WORK" --alpha "$A" --params "$GEOM/case_params.json"

cd "$WORK" || exit 1
simpleFoam > log.simpleFoam 2>&1
tail -1 postProcessing/forceCoeffs1/0/coefficient.dat > "coeff.txt" 2>/dev/null
echo "$GID a=$A node=$(hostname) $(cat coeff.txt 2>/dev/null)"
