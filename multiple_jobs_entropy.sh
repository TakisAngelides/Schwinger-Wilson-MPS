#!/bin/bash
echo "RUNNING $0"
N=$1
if [ -z $N ]; then
echo "MISSING VALUE OF N"
exit 1
fi
mg_list=(-0.3)
x=$2
if [ -z $x ]; then
echo "MISSING VALUE OF x"
exit 1
fi

for mg in "${mg_list[@]}"; do
	sbatch -J TN_${N}_${x}_${mg}_${D} job_Schwinger_entropy.sh ${N} ${x} ${mg}
done