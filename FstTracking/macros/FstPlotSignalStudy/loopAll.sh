#!/bin/bash

#. ./loopAll.sh

if [ $# -eq 0 ]
then
  mod=Mod03

  InPutHV=../../list/FST/data/${mod}/HighVoltage_${mod}.list
  for hv in `cat $InPutHV`
  do
    echo "Plot Noise Study for ${hv}"
    . ./plotQA.sh ${mod} ${hv}
  done
fi
