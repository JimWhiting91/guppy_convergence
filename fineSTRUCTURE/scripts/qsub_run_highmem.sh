#!/bin/bash


function show_help {
    echo "Usage: ./qsub_run.sh -f <cmdfile>"
    echo "  Submit commands to a qsub system, where <cmdfile> is a file containing a list of commands to execute, one per line."
    echo "OPTIONS:"
    echo "  -f <val>: The name of the command file to submit. REQUIRED."
    echo "  -n <val>: Set the number of processes per node requested in the qsub script. (Default: 1)"
    echo "  -w <val>: Set the walltime, in hours, requested in the qsub script. (Default: 120)"
    echo "  -s <val>: Set the sleep time in seconds between job submissions, which prevents the submission system from being overwhelmed. (Default: 10)"
    echo "  -m <val>: Set the maximum number of commands per qsub job. Using a large value is significantly more efficient than running each individually, provided you stay under the walltime. (Default: number of nodes)"
    echo "  -q <val>: Change the command run to process the script. This can be used for some non-qsub systems. (Default: \"qsub\")"
    echo "  -p      : Pretend: do not actually run the commands, but generate the qsub scripts for examination."
    echo "  -P <val>: Parallel execution tool. If set, all commands in a job will be written to a file and called via <val>. \"-P parallel\" is highly recommended if you have installed GNU parallel. If you have not, you can still use this tool with -P \"\" but -m will not work, and you MUST end each command with & for it to run in parallel! See http://www.gnu.org/software/parallel .  Default: parallel" 
    echo "  -v      : verbose mode."
    echo "  -h      : This help."
    exit 0
}

function makefile {
    pbsfile="$dir/t_${cmdfile}_cmd${filenum}.pbs"
    echo "#!/bin/bash
#PBS -l walltime=$walltime:00:00
#PBS -l nodes=1:ppn=$nodes
#PBS -o $dir/t_${cmdfile}_cmd${filenum}.log
#PBS -e $dir/t_${cmdfile}_cmd${filenum}.err
#PBS -q pq 
#PBS -A Research_Project-T110748
#PBS -d .
#PBS -l feature=highmem

# I have added this so that the libgsl.so.0 symlink is in the library path. Links to ~./linuxbrew/libs/libgsl.so
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/software/fs_4.0.1

# Load modules
module load GLib/2.53.5-GCCcore-6.4.0
module load GSL/2.4-GCCcore-6.4.0
module load Perl/5.28.0-GCCcore-7.3.0
" > $curdir/$pbsfile
 #  chmod +x $curdir/$pbsfile
 # echo cd $curdir >> $curdir/$pbsfile

}

function makepfile {
    pbsfile="$dir/t_${cmdfile}_cmd${filenum}.pbs"
    parallelfile="$dir/t_${cmdfile}_cmd${filenum}.parallel"
    rm -f $parallelfile
    touch $parallelfile
    echo "cd ../../" >> $pbsfile
    echo "cat $curdir/$parallelfile | $parallel -j $nodes" >> $pbsfile
}

function appendpfile {
    parallelfile="$dir/t_${cmdfile}_cmd${filenum}.parallel"
    echo "$@" >> $parallelfile
}

cmdfile=""
qsub="msub"
cmdsperproc=-1
verbose=0
sleepqsub=10
walltime=120  #hours
nodes=0
parallel="parallel"
pretend=0
while getopts "h?vpP:n:q:s:w:m:f:" opt; do
    case "$opt" in
    f) cmdfile=$OPTARG
	;;
    m)  cmdsperproc=$OPTARG
        ;;
    q)  qsub=$OPTARG
        ;;
    s) sleepqsub=$OPTARG
        ;;
    w) walltime=$OPTARG
        ;;
    n) nodes=$OPTARG
        ;;
    P) parallel=$OPTARG
        ;;
    p)  pretend=1
        ;;
    h|\?)
        show_help
        ;;
    v)  verbose=1
        ;;
    esac
done

if [ "$parallel" == "" ]; then
    echo "WARNING: Running without parallel. Commands will be exectuted as-is, as if in a bash file. They will ONLY run in parallel if run in the background with & at the end of each file. See qsub_run.sh -h for download instructions."
else
    test=`which $parallel`
    if [ "$test" == "" ] ; then
	echo "ERROR: You are trying to use $parallel for parallelization but this could not be found. Either download GNU parallel or see the -P option in the help."
	show_help
	exit 1
    fi
fi

if [ "$nodes" -eq "0" ]; then
    echo "ERROR: You did not set the number of processes per node. Set it with -n <n>. Recommend using 8, or 1/2 of the number of cores in the majority of nodes." 
    show_help
    exit 1
fi

if [ "$cmdfile" == "" ]; then
    echo "Command file is not set!"
    show_help
    exit 1
fi

if [ ! -e "$cmdfile" ]; then
    show_help
fi
if [[ $cmdsperproc -le 0 ]] ; then
    cmdsperproc=$nodes
fi

ncmds=`wc -l $cmdfile | cut -f1 -d' '`
if [ $ncmds -eq 0 ]; then
    echo "Command file $cmdfile appears to be empty!"
    exit 1
fi

curdir=`pwd`
ext=${cmdfile##*.}
dir=`basename $cmdfile .$ext`
dir="tmp_$dir"
mkdir -p "$dir"

linenum=0
filenum=1
empty=1
makefile
if [[ $parallel != "" ]]; then makepfile; fi

while read cmd
do
    empty=0
    if [[ $parallel != "" ]]; then 
	appendpfile "$cmd"; 
    else
	echo "$cmd" >> $curdir/$pbsfile
    fi
    if [ "$verbose" -eq 1 ]; then
	echo "   $pbsfile: $cmd"
    fi
    linenum=$(( $linenum + 1 ))
# if we've reached enough commands to finish
    if [ `echo $linenum $cmdsperproc |awk '{print ($1 % $2 == 0)}'` -eq 1 ] ;then
## submit the job
	echo $qsub $pbsfile
	if [ "$pretend" -eq 0 ]; then
	    $qsub $pbsfile
	fi
## Create an empty new job
	empty=1
	filenum=$(( $filenum + 1 ))
	makefile
	if [[ $parallel != "" ]]; then makepfile; fi
        sleep $sleepqsub
    fi
done < "$cmdfile"
if [ $empty -eq 0 ]; then
    echo $qsub $pbsfile
    if [ "$pretend" -eq 0 ]; then
	$qsub $pbsfile
    fi
fi
