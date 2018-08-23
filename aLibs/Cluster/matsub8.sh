# rm scratch/*
export NUM_NODES=8
qsub -N 0e -M fjs2@yale.edu -j oe -q scavenge -l nodes=1:ppn=16 -t 0-7 -v NUM_NODES re3run.sh

