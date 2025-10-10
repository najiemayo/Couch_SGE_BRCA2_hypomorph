#!/bin/sh


codeBase="$1"
outputRoot="$2"
params="$3"
funR="$4"
runsetf="$5"

cd ${codeBase}

source $runsetf

queue="cpu-short" ##"cpu-long""cpu-short"

mkdir -p ${outputRoot}/jobStats/

for exons in "${exonsets[@]}"; do
  for loessSubset in "${loessSubsets[@]}"; do
    for EventCount_Filter in ${EventCount_Filters[@]}; do
      for Ratio in ${Ratios[@]}; do
        for vset in ${vsets[@]}; do
          for fc in ${fcs[@]}; do
            
          parameterID=${exons}_${loessSubset}_${EventCount_Filter}_${Ratio}_${vset}_${fc}
          echo $parameterID
          rmdnm=${outputRoot}/jobStats/$parameterID.Rmd
          cp -f ./norm_classification_parm.Rmd $rmdnm
          sleep 2
          
          sbatch -p $queue -t 4:00:00 --mem=50G \
          -o ${outputRoot}/jobStats/$parameterID.out \
          -e ${outputRoot}/jobStats/$parameterID.err \
           ./runXReport.R $exons $loessSubset $EventCount_Filter $Ratio $vset $fc $rmdnm $params $funR
          done
        done
      done
    done
  done
done