# rm scratch/*
export NUM_JOBS=16
qsub -N 0e -M fjs2@yale.edu -j oe -q scavenge -l nodes=1:ppn=8 -t 0-7 -W PARTITION:m610 -v NUM_JOBS re3run2.sh

