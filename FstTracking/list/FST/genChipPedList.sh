#!/bin/bash
date

#. ./genChipPedList.sh Mod03 Inner/Outer

if [ $# -eq 2 ]
then
  mod=$1
  sector=$2

  mkdir -p noise/${mod}

  rm noise/${mod}/FstChipPed_${mod}_${sector}.list
  touch noise/${mod}/FstChipPed_${mod}_${sector}.list
  ls ../../data/FstCosmicTestStand_*/ped/Fst${sector}Chip*${mod}*_ped_*.root >> noise/${mod}/FstChipPed_${mod}_${sector}.list
fi
