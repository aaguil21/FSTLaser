#!/bin/bash
date

#. ./fstCluster.sh Mod03 hv isApplyCMNCorr nFstHitsCut date

if [ $# -eq 5 ]
then
  mod=$1
  hv=$2
  isSavePed=true
  isApplyCMNCorr=$3
  nFstHitsCut=$4
  date=$5

  root -l -b -q makeFstCluster.C\(\"${mod}\",\"${hv}\",${isSavePed},${isApplyCMNCorr},${nFstHitsCut},2,2.5,3.5,\"${date}\"\)
fi
