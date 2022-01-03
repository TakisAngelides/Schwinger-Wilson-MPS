#!/bin/bash
#SBATCH -o /onyx/qdata/TakisAngelides/logs/%x_%j.out
#SBATCH -e /onyx/qdata/TakisAngelides/logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --time=24:00:00
#SBATCH --partition=phi
#SBATCH --mem=31GB

N=$1
D=$2
D_previous=$3
x=$4
mg=$5

#/onyx/qdata/julia-1.7.0/bin/julia run_Schwinger_Wilson.jl 100 100 80 80.0 0.125
/onyx/qdata/julia-1.7.0/bin/julia run_Schwinger_Wilson.jl $1 $2 $3 $4 $5
#/onyx/qdata/julia-1.7.0/bin/julia testing.jl
