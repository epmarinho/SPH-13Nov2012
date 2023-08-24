#!/bin/bash
simcode=$1
seriesname=$2
m=$3
n=$4
N=$5
K=$6
omega=$7
epsilon=$8
dt=$9
if [ -e $seriesname$m.data.gz ]
then
	gunzip $seriesname$m.data.gz
	gzflag=1
fi
for((i=m;i<n;i++))
do
	let j=i+1
	time $simcode $seriesname$i.data $N $K $omega $epsilon $dt >$seriesname$j.data
	echo $i '->' $j
	[ -z $gzflag ] || gzip $seriesname$i.data
done
[ -z $gzflag ] || gzip $seriesname$m.data
