#!/bin/bash
date

#. ./genAlignmentList.sh Mod03

if [ $# -eq 1 ]
then
  mod=$1
  nFstHitsCut=4.0
  numOfUsedTimeBins=2
  config=withPed_withCMNCorr

  mkdir -p ./alignment/${mod}

  OutPutList=./alignment/${mod}/FstAlignment_${mod}.list
  rm $OutPutList
  touch $OutPutList

  InPutHV=data/${mod}/HighVoltage_${mod}.list
  for hv in `cat $InPutHV`
  do
    echo "Put ${hv} into Alignment List"
    InPutDate=../../list/FST/data/${mod}/Date_${mod}_${hv}.list
    InPutList=cluster/${mod}/FstCluster_${mod}_${hv}.list
    for item in `cat $InPutDate`
    do
      cat $InPutList | grep Th${nFstHitsCut}Tb${numOfUsedTimeBins}Ped2.5Ped3.5_${config}_${item}
      cat $InPutList | grep Th${nFstHitsCut}Tb${numOfUsedTimeBins}Ped2.5Ped3.5_${config}_${item} >> $OutPutList
    done
  done
fi
