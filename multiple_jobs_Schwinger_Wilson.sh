#!/bin/bash
echo "RUNNING $0"
D=$1
if [ -z $D ]; then
echo "MISSING VALUE OF D"
exit 1
fi
<<<<<<< HEAD
N_val=(40)
x_val=(30.0 40.0 50.0 60.0 80.0)
#N_val=(40)
#x_val=(30.0)
=======
N_val=(40 60 80 100)
x_val=(30.0 40.0 50.0 60.0 80.0)
>>>>>>> 018c0c40bdd86a45803666810df4613c5d735d57
mg=0.125
D_previous=$2
if [ -z $D_previous ]; then
echo "MISSING VALUE OF D_previous"
exit 1
fi


for N in "${N_val[@]}"; do 
	for x in "${x_val[@]}"; do 
<<<<<<< HEAD
		sbatch -J TN_${N}_${D}_${D_previous}_${x}_${mg} job_Schwinger_Wilson.sh ${N} ${D} ${D_previous} ${x} ${mg}
=======
		echo sbatch -J TN_${N}_${D}_${D_previous}_${x}_${mg} job_Schwinger_Wilson.sh ${N} ${D} ${D_previous} ${x} ${mg}
>>>>>>> 018c0c40bdd86a45803666810df4613c5d735d57
done
done
