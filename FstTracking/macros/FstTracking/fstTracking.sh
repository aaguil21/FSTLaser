#!/bin/bash
date

#. ./fstTracking.sh mod isApplyCMNCorr nFstHitsCut

if [ $# -eq 2 ]
then
  mod=$1
  hv=$2
  isSavePed=true
  isApplyCMNCorr=true
  nFstHitsCut=4.0

  root -l -b -q doFstTracking.C\(\"${mod}\",\"${hv}\",${isSavePed},${isApplyCMNCorr},${nFstHitsCut},2,2.5,3.5\)
fi
