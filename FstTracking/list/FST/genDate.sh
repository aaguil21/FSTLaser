#!/bin/bash
date

#. ./genDate.sh Mod03 HV200V

if [ $# -eq 2 ]
then
  mod=$1
  hv=$2

  rm data/${mod}/Date_${mod}_${hv}.list
  touch data/${mod}/Date_${mod}_${hv}.list
  cut -d '_' -f 5 data/${mod}/FstData_${mod}_${hv}.list | sort | uniq >> data/${mod}/Date_${mod}_${hv}.list
fi
