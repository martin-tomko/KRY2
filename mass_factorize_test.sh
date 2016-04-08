#!/bin/bash

gen="./semiprimegen.py"
timeouts="cat timeouts.txt"

prog=./kry
prog_temp=$prog-$(date +%H%M%S)

cmd=$gen
time_cmd="/usr/bin/time -f %U"
time_limit=10
verify_cmd="./verify_factor.py"

cp $prog $prog_temp

for num in `$cmd`; do
  timeout $time_limit $time_cmd $prog_temp -b $num || echo "Timeout at $num" >&2
#  result=`$prog_temp -b $num`
#  $verify_cmd $num $result || echo "Bad: $num $result" >&2
done

rm $prog_temp
