# matsub.sh nNodes PPN nodeType script
# rm scratch/*
if [ $# -lt 3 ]
then
	echo usage:
	echo "matsub.sh <nNodes> <PPN> <nodeType> <script>"
	echo e.g.
	echo matsub.sh 4 16 m620 matlocalk2p.sh
else
	rm 0i*
	export NUM_NODES=$1
	qsub -N 0i -M fjs2@yale.edu -j oe -q scavenge -l nodes=1:ppn=$2 -t 0-$((NUM_NODES - 1)) -W PARTITION:$3 -v NUM_NODES $4
fi

