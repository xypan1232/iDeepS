#!/bin/bash
#$ -cwd
#$ -l h_vmem=2G
#$ -sync yes

######################################################################
## job script "generic_submit_to_cluster.sh"
## so this is a generic shell script that can call a certain perl script
## with the given parameters. 
## This is useful if you need to commit several similar jobs to the cluster
######################################################################

#standard path variable for any programs that are used
export PATH=/usr/local/vrna/1.8.5/bin/:$HOME/bin/:$PATH

# all parameters for the shell script
working_dir=$1
perlscript=$2
params=$3

echo working_dir: $working_dir
echo perlscript: $perlscript
echo params: $params

if [ -z "$SGE_TASK_ID" ]; then
	echo "No SGE environment found! Set SGE_TASK_ID to 0!"
	let SGE_TASK_ID=0
fi

echo call: "/usr/local/perl/bin/perl $perlscript $params --jobid $SGE_TASK_ID"

cd $working_dir
/usr/local/perl/bin/perl $perlscript $params --jobid $SGE_TASK_ID

