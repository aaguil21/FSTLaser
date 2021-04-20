#!/bin/bash

#. ./loopAll.sh

if [ $# -eq 0 ]
then
  isApplyCMNCorr=false

  echo "Time Sequence Study for HV140V"
  InPutDate_HV140V=../../list/FST/DateHV140V.list
  for item in `cat $InPutDate_HV140V`
  do
   . ./fstQAStudy_timesequence.sh HV140V $isApplyCMNCorr $item 
  done

  echo "Time Sequence Study for HV200V"
  InPutDate_HV200V=../../list/FST/DateHV200V.list
  for item in `cat $InPutDate_HV200V`
  do
   . ./fstQAStudy_timesequence.sh HV200V $isApplyCMNCorr $item 
  done
fi
