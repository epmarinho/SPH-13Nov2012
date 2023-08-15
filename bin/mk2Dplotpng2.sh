#!/bin/sh
objname="$1"
for f in $objname*.d
do
	(echo set term 'png' size 512,512;
	 echo set output \"`basename $f d`png\";
	 echo plot[0:1.1][-25:0]\"$f\"u \(\$2\*\$2+\$3\*\$3+\$4\*\$4\):\(log\(\$8\)\) w d notitle
	)|gnuplot
done
