#!/bin/bash
date

#. ./splitCluster_allConfig.sh hv

if [ $# -eq 1 ]
then
  hv=$1

  echo "Split Cluster list for ${hv}"  
   . ./splitCluster_config.sh ${hv} 4.5 3 withPed_withCMNCorr
   . ./splitCluster_config.sh ${hv} 4.5 2 withPed_withCMNCorr
   # . ./splitCluster_config.sh ${hv} 4.5 1 withPed_withCMNCorr
   . ./splitCluster_config.sh ${hv} 4.0 3 withPed_withCMNCorr
   . ./splitCluster_config.sh ${hv} 4.0 2 withPed_withCMNCorr
   # . ./splitCluster_config.sh ${hv} 4.0 1 withPed_withCMNCorr
   # . ./splitCluster_config.sh ${hv} 3.5 3 withPed_withCMNCorr
   # . ./splitCluster_config.sh ${hv} 3.5 2 withPed_withCMNCorr
   # . ./splitCluster_config.sh ${hv} 3.5 1 withPed_withCMNCorr
   # . ./splitCluster_config.sh ${hv} 3.0 3 withPed_withCMNCorr
   # . ./splitCluster_config.sh ${hv} 3.0 2 withPed_withCMNCorr
   # . ./splitCluster_config.sh ${hv} 3.0 1 withPed_withCMNCorr

   . ./splitCluster_config.sh ${hv} 4.5 3 withPed_woCMNCorr
   . ./splitCluster_config.sh ${hv} 4.5 2 withPed_woCMNCorr
   # . ./splitCluster_config.sh ${hv} 4.5 1 withPed_woCMNCorr
   . ./splitCluster_config.sh ${hv} 4.0 3 withPed_woCMNCorr
   . ./splitCluster_config.sh ${hv} 4.0 2 withPed_woCMNCorr
   # . ./splitCluster_config.sh ${hv} 4.0 1 withPed_woCMNCorr
   # . ./splitCluster_config.sh ${hv} 3.5 3 withPed_woCMNCorr
   # . ./splitCluster_config.sh ${hv} 3.5 2 withPed_woCMNCorr
   # . ./splitCluster_config.sh ${hv} 3.5 1 withPed_woCMNCorr
   # . ./splitCluster_config.sh ${hv} 3.0 3 withPed_woCMNCorr
   # . ./splitCluster_config.sh ${hv} 3.0 2 withPed_woCMNCorr
   # . ./splitCluster_config.sh ${hv} 3.0 1 withPed_woCMNCorr
fi
