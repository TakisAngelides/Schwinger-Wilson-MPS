#!/bin/bash
echo "RUNNING $0"
N=$1
x=$2
D=$3
if [ -z $N ]; then
echo "MISSING VALUE OF N"
exit 1
fi
if [ -z $x ]; then
echo "MISSING VALUE OF x"
exit 1
fi
if [ -z $D ]; then
echo "MISSING VALUE OF x"
exit 1
fi
mg_list=(-0.1 -0.2 -0.3 -0.45 -0.47 -0.49 -0.51 -0.53 -0.55 -0.6 -0.65 -0.75)

for mg in "${mg_list[@]}"; do
	sbatch -J TN_${N}_${x}_${mg}_${D} job_Schwinger_entropy.sh ${N} ${x} ${mg} ${D}
done