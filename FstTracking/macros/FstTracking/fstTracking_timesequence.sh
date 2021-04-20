#!/bin/bash
date

#. ./fstTracking_timesequence.sh hv isApplyCMNCorr date

if [ $# -eq 3 ]
then
  hv=$1
  isSavePed=true
  isApplyCMNCorr=$2
  date=$3

  root -l -b -q doFstTracking_timesequence.C\(\"${hv}\",${isSavePed},${isApplyCMNCorr},4.5,3,2.0,\"${date}\"\)
fi
