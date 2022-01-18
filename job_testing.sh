#!/bin/bash
#SBATCH -o /onyx/qdata/TakisAngelides/logs/%x_%j.out
#SBATCH -e /onyx/qdata/TakisAngelides/logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --time=24:00:00
#SBATCH --partition=phi
#SBATCH --mem=25GB

/onyx/qdata/julia-1.7.0/bin/julia testing.jl