#!/bin/bash

#. ./loopAll.sh

if [ $# -eq 0 ]
then
  mod=Mod04
  isApplyCMNCorr=true
  # isApplyCMNCorr=false

  InPutHV=../../list/FST/data/${mod}/HighVoltage_${mod}.list
  for hv in `cat $InPutHV`
  do
    echo "Find Clusters for ${hv}"
    InPutDate=../../list/FST/data/${mod}/Date_${mod}_${hv}.list
    for item in `cat $InPutDate`
    do
      . ./fstCluster.sh ${mod} ${hv} $isApplyCMNCorr 4.0 $item 
      # . ./fstCluster.sh ${mod} ${hv} $isApplyCMNCorr 3.5 $item 
    done
  done
fi
