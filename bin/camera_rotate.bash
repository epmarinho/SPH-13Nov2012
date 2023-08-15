#!/bin/bash
basename=$1
thetax=$2
thetay=$3
thetaz=$4
delta=$5
m=$6
n=$7
rotate=/home/projects/spmd_treecode/src/rotate
let n++
for((i=m;i<n;i++))
do
  dfile=$basename$i.data
  echo $dfile
  gunzip -f $dfile.gz
  $rotate $thetax $thetay $thetaz < $dfile > rotate-$dfile
  thetay=`echo $thetay + $delta | /usr/bin/bc -l`
  gzip -f $dfile
  gzip -f rotate-$dfile
done
