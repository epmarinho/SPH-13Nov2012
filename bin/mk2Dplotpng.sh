#!/bin/sh
objname="$1"
a="$2"
b="$3"
c="$4"
d="$5"
x="$6"
y="$7"
sx="$8"
sy="$9"
for f in $objname*.d
do
	(echo set term 'png' size $sx,$sy;
	 echo set output \"`basename $f d`png\";
	 echo plot[$a:$b][$c:$d]\"$f\"u $x:$y w d notitle
	)|gnuplot
done
