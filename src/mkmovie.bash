#!/bin/bash
#@ mkmovie.bash::
#@ Eraldo Pereira Marinho
#@ Mon Sep 17 18:54:53 BRT 2012
#@ Version 1
#@ Revision 1.2
#@ Revision 1.3
objname="$1"
a="$2"
b="$3"
c="$4"
d="$5"
x="$6"
y="$7"
sx="$8"
Sx="${9}"
m="${10}"
n="${11}"
dcolor="${12}"
plotxy=~/projects/sph/bin/plotxy$dcolor
let n++
for((i=m;i<n;i++))
do
	let j=100000+i
	datafile=$objname$i.data
	pngfile=$objname$j.png
	if [ ! -e $pngfile ]
	then
		echo $datafile '->' $pngfile
		if [ -e $datafile.gz ]
		then
			gunzip $datafile.gz
			gzflag=1
		else
			gzflag=
		fi
		rasterfile=$objname$j.ras
		$plotxy $datafile $a $b $c $d $x $y $sx > $rasterfile
		[ -z $gzflag ] || gzip -f $datafile
		convert -resize ${Sx}x$Sx $rasterfile $pngfile
		rm -f $rasterfile
	fi
done
mencoder "mf://$objname*.png" -noskip -mf fps=60 -o $objname.avi -ovc lavc -lavcopts vcodec=mpeg4
