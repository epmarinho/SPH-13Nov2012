#!/bin/sh
objname="$1"
s=$2
for f in $objname[1-9]*.d
do
	(echo set term 'png' size 400,550;
	 echo set output \"`basename $f d`png\";
	 echo splot[-$s:$s][-$s:$s][-$s:$s]\"planet.d\"u 2:3:4 w d lc rgb \'orange\',\"$f\"u 2:3:4 w d lc rgb \'blue\' notitle
	)|gnuplot
done
