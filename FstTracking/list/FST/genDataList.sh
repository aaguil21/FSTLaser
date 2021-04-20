#!/bin/bash
date

#. ./genDataList.sh Mod03 HV200V

if [ $# -eq 2 ]
then
  mod=$1
  hv=$2

  mkdir -p data/${mod}

  rm data/${mod}/FstData_${mod}_${hv}.list
  touch data/${mod}/FstData_${mod}_${hv}.list
  ls ../../data/FstCosmicTestStand_*/data/Fst*${mod}*_${hv}_data_*.root >> data/${mod}/FstData_${mod}_${hv}.list
fi
