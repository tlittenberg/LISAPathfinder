#!/bin/sh

#This is to propagate interrupt to child processes
trap 'killall' INT
killall() {
    trap '' INT TERM     # ignore INT and TERM while shutting down
    echo "**** Shutting down... ****"     # added double quotes
    kill -TERM 0         # fixed order, send TERM not INT
    wait
    echo DONE
}


ncount=$1
nproc=$2
BASEDIR=$(dirname "$0")

for i in `seq 1 $nproc`;
do
    $BASEDIR/pprun1.sh $ncount $i &
done
wait
