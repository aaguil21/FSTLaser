#!/bin/bash
# date

#. ./splitCluster.sh Mod03 HV70V
if [ $# -eq 2 ]
then
  mod=$1
  hv=$2
  nFstHitsCut=4.0
  numOfUsedTimeBins=2
  config=withPed_withCMNCorr

  InPutList=cluster/${mod}/FstCluster_${mod}_${hv}.list
  InPutDate=data/${mod}/Date_${mod}_${hv}.list

  OutPutList=./cluster/${mod}/FstCluster_${mod}_${hv}_Th${nFstHitsCut}Tb${numOfUsedTimeBins}Ped2.5Ped3.5_${config}.list
  echo "OutPut list is: ${OutPutList}"
  rm $OutPutList
  touch $OutPutList
  for item in `cat $InPutDate`
  do
    # echo $item
    cat $InPutList | grep Th${nFstHitsCut}Tb${numOfUsedTimeBins}Ped2.5Ped3.5_${config}_${item}
    cat $InPutList | grep Th${nFstHitsCut}Tb${numOfUsedTimeBins}Ped2.5Ped3.5_${config}_${item} >> $OutPutList
  done
fi
