#!/bin/bash
date

#. ./genPedHighVoltage.sh Mod01

if [ $# -eq 1 ]
then
  mod=$1

  rm noise/${mod}/HighVoltage_${mod}.list
  touch noise/${mod}/HighVoltage_${mod}.list
  cut -d '_' -f 6 noise/${mod}/FstPed_${mod}_HV*.list | sort | uniq >> noise/${mod}/HighVoltage_${mod}.list
fi
