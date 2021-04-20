#!/bin/bash
date

#. ./fstQAStudy.sh Mod03 HV70V

if [ $# -eq 2 ]
then
  mod=$1
  hv=$2
  isSavePed=true
  isApplyCMNCorr=true
  nFstHitsCut=4.0

  root -l -b -q doFstQAStudy.C\(\"${mod}\",\"${hv}\",${isSavePed},${isApplyCMNCorr},${nFstHitsCut},2,2.5,3.5\)
fi
