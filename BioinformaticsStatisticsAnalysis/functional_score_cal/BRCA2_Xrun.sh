#!/bin/sh

## before run this chmod +x this file

codeBase="rpath/MAVE/"
outputRoot="rpath/pipeline/"
params="/rpath/param_BRCA2.R"
funR="/rpath/MAVE/functions.R"
runsetf="/rpath/BRCA2_Xrun_params.txt"

${codeBase}/SGErun.sh $codeBase $outputRoot $params $funR $runsetf
