#!/bin/sh

ncount=$1
pnum=$2
BASEDIR=$(dirname "$0")
EXEC=${BASEDIR}/../mcmc
base=lock
echo "entering pprun1 on proc $pnum: ncount=$ncount"
for i in `seq 1 $ncount`;
do
    file="$base$i"
    if ( set -o noclobber; echo "$pnum" > "$file") 2> /dev/null; then
	trap 'rm -f "$lockfile"; exit $?' INT TERM EXIT
	echo "Locking succeeded" >&2
	echo "proc $pnum doing $i"
	rundir="run$i"
	mkdir -p $rundir
	cd $rundir
	echo ../${EXEC} -d 3 -s $RANDOM -n $RANDOM -j >mcmc.out
	../${EXEC} -d 3 -s $RANDOM -n $RANDOM -j >> mcmc.out 2>&1
	#sleep `bc -l <<< "$RANDOM / 262144 "`
	cd ..
    fi
done

