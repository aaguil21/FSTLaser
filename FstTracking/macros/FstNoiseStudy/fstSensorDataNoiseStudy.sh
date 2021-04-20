#!/bin/bash
date

#. ./fstSensorDataNoiseStudy.sh

if [ $# -eq 3 ]
then
  mod=$1
  hv=$2
  date=$3

  root -l -b -q doFstSensorDataNoiseStudy.C\(\"${mod}\",\"${hv}\",\"${date}\"\)
fi
