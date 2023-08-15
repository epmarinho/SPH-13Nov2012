#!/bin/bash
simcode=$1
seriesname=$2
m=$3
n=$4
omega=$5
epsilon=$6
dt=$7
if [ -e $seriesname$m.data.gz ]
then
	gunzip $seriesname$m.data.gz
	gzflag=1
fi
for((i=m;i<n;i++))
do
	let j=i+1
	time $simcode $seriesname$i.data $omega $epsilon $dt >$seriesname$j.data
	echo $i '->' $j
	[ -z $gzflag ] || gzip $seriesname$i.data
done
[ -z $gzflag ] || gzip $seriesname$m.data
