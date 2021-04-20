#!/bin/bash
date

#. ./plotSignalTS.sh

if [ $# -eq 0 ]
then
  # root -l -b -q plotSignalTS_HV70V.C\(\)
  # root -l -b -q plotSignalTS_HV70V_APV.C\(0,0\)
  # root -l -b -q plotSignalTS_HV70V_APV.C\(0,1\)
  # root -l -b -q plotSignalTS_HV70V_APV.C\(1,0\)
  # root -l -b -q plotSignalTS_HV70V_APV.C\(1,1\)

  root -l -b -q plotSignalTS_HV140V.C\(\)
  root -l -b -q plotSignalTS_HV140V_APV.C\(0,0\)
  root -l -b -q plotSignalTS_HV140V_APV.C\(0,1\)
  root -l -b -q plotSignalTS_HV140V_APV.C\(1,0\)
  root -l -b -q plotSignalTS_HV140V_APV.C\(1,1\)

  root -l -b -q plotSignalTS_HV200V.C\(\)
  root -l -b -q plotSignalTS_HV200V_APV.C\(0,0\)
  root -l -b -q plotSignalTS_HV200V_APV.C\(0,1\)
  root -l -b -q plotSignalTS_HV200V_APV.C\(1,0\)
  root -l -b -q plotSignalTS_HV200V_APV.C\(1,1\)
fi
