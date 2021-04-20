#!/bin/bash

#. ./loopAll.sh

if [ $# -eq 0 ]
then
  mod=Mod04

  InPutHV=../../list/FST/data/${mod}/HighVoltage_${mod}.list
  for hv in `cat $InPutHV`
  do
    echo "Do QA Study for ${hv}"
    . ./fstQAStudy.sh ${mod} ${hv}
  done
fi
