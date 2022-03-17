#!/bin/sh

N=$1
x=$2
D=$3

for mg in {-0.1 -0.2 -0.3 -0.45 -0.47 -0.49 -0.51 -0.53 -0.55 -0.6 -0.65 -0.75}
do
    /lustre/fs23/group/nic/tangelides/julia-1.7.2/bin/julia /afs/ifh.de/user/a/angeltak/Schwinger-Wilson-MPS/run_Schwinger_Wilson.jl $N $x $mg $D
done