#!/bin/bash

#. ./loopSensorDataNoise.sh

if [ $# -eq 0 ]
then
  mod=Mod03

  InPutHV=../../list/FST/noise/${mod}/HighVoltage_${mod}.list
  for hv in `cat $InPutHV`
  do
    echo "Do Noise Scan for ${hv}"
    . ./fstSensorPedNoiseStudy.sh ${mod} ${hv}
  done
fi
