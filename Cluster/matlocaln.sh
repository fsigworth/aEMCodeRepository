# matlocaln numNodeJobs scriptName
#  launch numNodeJobs instances of a script.
#  defaults for arguments are numNodeJobs=8, scriptName = reEMReconstruct3
# We pass to Matlab the following environment variables:
# BATCH_ENVIRONMENT = BATCH
# NUM_NODES (hopefully already set)
# NUM_JOBS = numNodeJobs * NUM_NODES
# PBS_ARRAY_ID (this is the job number, runs from 0...
#                      numNodeJobs * originalPBS_ARRAY_ID-1
# So if we launch a PBS array, this is expanded into these intance numbers.
# The lowest PBS_ARRAY_ID value is passed to the blocking job.
#
echo matlocaln.sh starting matlab jobs.

matlab="/Applications/MATLAB_R2015b.app/bin/matlab -nojvm"

if [ $# -gt 0 ]
then
	num=$1
else
	num=4 # default number of node jobs
fi

if [ $# -gt 1 ]
then
	exec=$2
else
	exec=reEMReconstruct3 # default execution script
fi
echo script = $exec
if [ $(( NUM_NODES )) -eq 0 ]
then
	NUM_NODES=1
fi
export NUM_JOBS=$((num * NUM_NODES))
export BATCH_ENVIRONMENT=BATCH
id0=$PBS_ARRAYID
i=0
while [ $((++i)) -lt $num ]
do
	export PBS_ARRAYID=$((i + num * id0))
	$matlab -r $exec &
	echo Job $PBS_ARRAYID
done
 First job is blocking
export PBS_ARRAYID=$((num * id0))
$matlab -r $exec
echo Job $PBS_ARRAYID

