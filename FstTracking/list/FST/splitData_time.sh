#!/bin/bash
date

#. ./splitData_time.sh Mod03 HV200V
if [ $# -eq 2 ]
then
  mod=$1
  hv=$2
  InPutList=data/${mod}/FstData_${mod}_${hv}.list
  InPutDate=data/${mod}/Date_${mod}_${hv}.list

  for item in `cat $InPutDate`
  do
    echo $item
    OutPutList=./data/${mod}/FstData_${mod}_${hv}_${item}.list
    rm $OutPutList
    touch $OutPutList
    cat $InPutList | grep $item >> $OutPutList
  done
fi
