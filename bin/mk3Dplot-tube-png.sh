#!/bin/sh
objname="$1"
for file in $objname[1-9]*.d
do
	(echo set term 'png' size 480,768
	 echo set output \"`basename $file d`png\";
	 echo set view 70,120;
	 echo splot[-.01:1.01][-.051:.051][]\"$file\"u 4:2:3:7:5:6 w vec lc rgb \"blue\" notitle
	)|gnuplot
done
