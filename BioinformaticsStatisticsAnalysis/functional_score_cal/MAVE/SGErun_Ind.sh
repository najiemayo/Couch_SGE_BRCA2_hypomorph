#!/bin/sh


codeBase="$1"
outputRoot="$2"
params="$3"
funR="$4"
runsetf="$5"

source $runsetf

queue="cpu-short"

mkdir -p ${outputRoot}/jobStats/
cd ${outputRoot}/jobStats/

  for loessSubset in "${loessSubsetFiters[@]}"; do
    for EventCount_Filter in ${EventCount_Filters[@]}; do
      for ladjt in ${ladjts[@]}; do
        for vset in ${vsets[@]}; do
            
          parameterID=${loessSubset}_${EventCount_Filter}_${ladjt}_${vset}
          echo $parameterID
          rmdnm=$parameterID.Rmd
          rmdpath=${outputRoot}/jobStats/
          cp -f ${codeBase}/Ind_Exon.Rmd $rmdnm
          sleep 2
          
          sbatch -p $queue -t 4:00:00 --mem=50G \
          -o ${outputRoot}/jobStats/$parameterID.out \
          -e ${outputRoot}/jobStats/$parameterID.err \
           ${codeBase}/runIndReport.R $loessSubset $EventCount_Filter $ladjt $vset $rmdnm $params $funR $rmdpath
        
        done
      done
    done
  done
