#!/bin/bash

#. ./loopSensorDataNoise.sh

if [ $# -eq 0 ]
then
  mod=Mod04

  InPutHV=../../list/FST/data/${mod}/HighVoltage_${mod}.list
  for hv in `cat $InPutHV`
  do
    echo "Noise Study for ${hv}"
    InPutDate=../../list/FST/data/${mod}/Date_${mod}_${hv}.list
    for item in `cat $InPutDate`
    do
      . ./fstSensorDataNoiseStudy.sh ${mod} ${hv} ${item}
    done
  done
fi
