#!/bin/bash

echo "$0: NOTE -- currently only runs very simple crypt/factor tests." >&2
echo "            Also, the factor results are only compared with the smallest prime." >&2

program="./kry"

function run_test () 
{ # arguments: command to run, expected result.
  cmd=$1
  expected=$2

#  echo "About to run: $cmd" >&2
#  echo "Expected result: $expected" >&2
#  echo "Actual result:   $($cmd)" >&2

  if [[ $($cmd) != $expected ]]; then
    echo 1
    echo "Failed at $cmd" >&2
    return 1
  else
    return 0
  fi
}

function test_crypt ()
{
  N=$1; E=$2; D=$3; M=$4; C=$5

  cmd="$program -e $E $N $M"
  run_test "$cmd" $C || return

  cmd="$program -d $D $N $C"
  run_test "$cmd" $M || return

  echo 0
}

function test_factor ()
{
  N=$1; P=$2

  cmd="$program -b $N"
  run_test "$cmd" $P || return

  echo 0
}





  # Initialization:

failed_cnt=0

  # Basic RSA encrypt/decrypt tests:

fail=`test_crypt 119 5 77 0x13 0x42`
failed_cnt=$(( failed_cnt + fail ))

fail=`test_crypt 77 37 13 0xf 0x47`
failed_cnt=$(( failed_cnt + fail ))


  # Factorization tests:

fail=`test_factor 8633 0x59`
failed_cnt=$(( failed_cnt + fail ))

fail=`test_factor 1689243484681 0x13d4fd`
failed_cnt=$(( failed_cnt + fail ))

fail=`test_factor 519376503391506929869 0x54dab50ad`
failed_cnt=$(( failed_cnt + fail ))

  # Final evaluation:

if [[ $failed_cnt -eq 0 ]]; then
  echo "$0: Passed all tests." >&2
  exit 0
else
  echo "$0: Failed $failed_cnt tests." >&2
  exit 1
fi

