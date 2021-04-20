#!/bin/bash
date

#. ./plotQA.sh Mod01 HV70V

if [ $# -eq 2 ]
then
  mod=$1
  hv=$2

  folder=./figures/${mod}
  mkdir -p ${folder} 

  # echo "Plot Ped Noise QA!!"
  root -l -b -q plotSensorPedNoiseQA.C\(\"${mod}\",\"${hv}\"\)

  # echo "Plot Data Noise QA!!"
  InPutDate=../../list/FST/data/${mod}/Date_${mod}_${hv}.list
  for item in `cat $InPutDate`
  do
    root -l -b -q plotSensorDataNoiseQA.C\(\"${mod}\",\"${hv}\",\"${item}\"\)
  done

  root -l -b -q plotSensorDataNoiseDate.C\(\"${mod}\",\"${hv}\"\)
fi
