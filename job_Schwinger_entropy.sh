#!/bin/bash
#SBATCH -o /onyx/qdata/TakisAngelides/logs/%x_%j.out
#SBATCH -e /onyx/qdata/TakisAngelides/logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --time=24:00:00
#SBATCH --partition=phi
#SBATCH --mem=31GB

# N=$1
# x=$2
# mg=$3

/onyx/qdata/julia-1.7.0/bin/julia run_Schwinger_Wilson.jl $1 $2 $3