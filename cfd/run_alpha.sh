#!/bin/sh
# One Slurm array task = one angle of attack. Builds a self-contained work
# dir that SYMLINKS the shared mesh (constant/polyMesh) so all 21 cases read
# one mesh, stamps the freestream/force directions, and runs simpleFoam serial.
BASE="$HOME/cfd/disc/cfd"
ALPHAS="-5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15"
IDX=$((${SLURM_ARRAY_TASK_ID:-0} + 1))
A=$(echo $ALPHAS | cut -d' ' -f $IDX)

WORK="$BASE/sweep/a_$A"
rm -rf "$WORK"; mkdir -p "$WORK/constant"
cp -r "$BASE/case/0" "$WORK/0"
cp -r "$BASE/case/system" "$WORK/system"
cp "$BASE/case/constant/transportProperties" "$BASE/case/constant/turbulenceProperties" "$WORK/constant/"
ln -s "$BASE/case/constant/polyMesh" "$WORK/constant/polyMesh"

cd "$BASE" || exit 1
python stamp_alpha.py --case "$WORK" --alpha "$A" --params case_params.json

cd "$WORK" || exit 1
simpleFoam > log.simpleFoam 2>&1
RC=$?
tail -1 postProcessing/forceCoeffs1/0/coefficient.dat > "RESULT_${A}.txt" 2>/dev/null
echo "alpha=$A rc=$RC node=$(hostname) $(cat RESULT_${A}.txt 2>/dev/null)"
