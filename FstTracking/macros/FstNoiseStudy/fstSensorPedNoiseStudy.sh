#!/bin/bash
date

#. ./fstSensorPedNoiseStudy.sh Mod01 HV70V

if [ $# -eq 2 ]
then
  mod=$1
  hv=$2

  root -l -b -q doFstSensorPedNoiseStudy.C\(\"${mod}\",\"${hv}\"\)
fi
