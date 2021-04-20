#!/bin/bash
date

#. ./genHighVoltage.sh Mod03

if [ $# -eq 1 ]
then
  mod=$1

  rm data/${mod}/HighVoltage_${mod}.list
  touch data/${mod}/HighVoltage_${mod}.list
  cut -d '_' -f 6 data/${mod}/FstData_${mod}_HV*V.list | sort | uniq >> data/${mod}/HighVoltage_${mod}.list
fi
