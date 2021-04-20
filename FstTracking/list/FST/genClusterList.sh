#!/bin/bash
date

#. ./genClusterList.sh Mod03 HV70V

if [ $# -eq 2 ]
then
  mod=$1
  hv=$2

  mkdir -p cluster/${mod}

  rm cluster/${mod}/FstCluster_${mod}_${hv}.list
  touch cluster/${mod}/FstCluster_${mod}_${hv}.list
  ls ../../output/clusters/${mod}/FstClusters_${mod}_${hv}_*.root >> cluster/${mod}/FstCluster_${mod}_${hv}.list
fi
