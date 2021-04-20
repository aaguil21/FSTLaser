#!/bin/bash
date

#. ./genPedList.sh Mod03 HV200V

if [ $# -eq 2 ]
then
  mod=$1
  hv=$2

  mkdir -p noise/${mod}

  rm noise/${mod}/FstPed_${mod}_${hv}.list
  touch noise/${mod}/FstPed_${mod}_${hv}.list
  ls ../../data/FstCosmicTestStand_*/ped/FstInOuterSensor*${mod}*_${hv}_ped_*.root >> noise/${mod}/FstPed_${mod}_${hv}.list
fi
