#!/bin/sh
# One Slurm array task = one (anchor disc, angle) pair. task_id = disc*7 + alpha_i.
BASE="$HOME/cfd/disc/cfd"; AN="$BASE/anchors"; TPL="$BASE/bigsweep/template"
ALPHAS="-4 -1 2 5 8 11 14"; NA=7
DISCS="anchor_aviar anchor_roc anchor_wraith"

T=${SLURM_ARRAY_TASK_ID:-0}
GI=$((T / NA)); AI=$((T % NA))
A=$(echo $ALPHAS | cut -d' ' -f $((AI + 1)))
GID=$(echo $DISCS | cut -d' ' -f $((GI + 1)))
GEOM="$AN/geoms/$GID"

[ -f "$GEOM/constant/polyMesh/owner" ] || { echo "no mesh for $GID (skip)"; exit 0; }

WORK="$AN/solve/${GID}_a${A}"
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
