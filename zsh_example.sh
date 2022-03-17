#!/bin/zsh
#
#(otherwise the default shell would be used) 
#$ -S /bin/bash
#(the running time for this job)
#$ -l h_rt=47:59:00
#
#(the maximum memory usage per core of this job)
#$ -l h_rss=2.0G
#
#(output for std_out)
#$ -o /lustre/fs23/group/nic/skuehn/spin-helix/Julia/data_more_statistics/delta-0.495/logs
#
#(output for std_err)
#$ -e /lustre/fs23/group/nic/skuehn/spin-helix/Julia/data_more_statistics/delta-0.495/logs
#
#(send mail on job's abort)
#$ -m a

# Export the julia path
export JULIA_DEPOT_PATH="/lustre/fs23/group/nic/skuehn/.julia:$JULIA_DEPOT_PATH"

# change to scratch directory
cd $TMPDIR

# create temporary directories
mkdir data
cd data

# Create the suffix from the index of the array job
suffix="_i${SGE_TASK_ID}"


# run the simulation
/lustre/fs23/group/nic/skuehn/julia-1.5.2/bin/julia /lustre/fs23/group/nic/skuehn/spin-helix/Julia/runNoisyEvoTobias.jl $1 $2 $suffix $3 $4 $5


# copy the output back into the lustre file system
cp fluctuating_evolution_*.txt /lustre/fs23/group/nic/skuehn/spin-helix/Julia/data_more_statistics2/delta-0.495/evodata
cp gamma_*.txt /lustre/fs23/group/nic/skuehn/spin-helix/Julia/data_more_statistics2/delta-0.495/gamma
cp spin_config_*.txt /lustre/fs23/group/nic/skuehn/spin-helix/Julia/data_more_statistics2/delta-0.495/config
cp spin_current_*.txt /lustre/fs23/group/nic/skuehn/spin-helix/Julia/data_more_statistics2/delta-0.495/current
cp statevectors_evolution_* /lustre/fs23/group/nic/skuehn/spin-helix/Julia/data_more_statistics2/delta-0.495/states