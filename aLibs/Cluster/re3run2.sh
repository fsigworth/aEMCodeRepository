id0=$PBS_ARRAYID
export PBS_ARRAYID=$((2*id0))
/usr/local/cluster/software/matlab-2014a/bin/matlab -r "reEMReconstruct3" &
export PBS_ARRAYID=$((2*id0+1))
/usr/local/cluster/software/matlab-2014a/bin/matlab -r "reEMReconstruct3"
