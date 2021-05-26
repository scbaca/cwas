#! /bin/awk -f
BEGIN {	s=0; }
{ s+=$1; }
END { print s; }
