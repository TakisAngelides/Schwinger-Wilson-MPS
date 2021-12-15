#!/bin/bash
#SBATCH -o /onyx/qdata/TakisAngelides/logs/%x_%j.out
#SBATCH -e /onyx/qdata/TakisAngelides/logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --time=23:00:00
#SBATCH --partition=phi
#SBATCH --mem-per-cpu=8G

/onyx/qdata/julia-1.7.0/bin/julia run_Schwinger_Wilson.jl 20 40