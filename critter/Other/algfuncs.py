# For CA-CQR2
camfs_cacqr2 () {
  # launch CQR2
  local scale=${1}
  local binaryPath=${2}
  local numIterations=${3}
  local launchID=${4}
  local NumNodes=${5}
  local ppn=${6}
  local tpr=${7}
  local matrixDimM=${8}
  local matrixDimN=${9}
  local matrixDimMorig=${10}
  local matrixDimNorig=${11}
  local pDimDorig=${12}
  local pDimCorig=${13}
  local pDimD=${14}
  local pDimC=${15}
  local nodeIndex=${16}
  local scaleRegime=${17}
  local nodeCount=${18}
  local WScounterOffset=${19}
  local tuneInvCutOff=${20}
  local bcDim=0
  local tag1="cacqr2"

  # Next: Based on pDimC, decide on invCutOff parameter, which will range from 0 to a max of 2 for now
  curInverseCutOffMult=0
  # Set up the file string that will store the local benchmarking results
  local fileString="DataFiles/results_${tag1}_${scale}_${NumNodes}nodes_${matrixDimM}dimM_${matrixDimN}dimN_${curInverseCutOffMult}inverseCutOffMult_0bcMult_0panelDimMult_${pDimD}pDimD_${pDimC}pDimC_${numIterations}numIter_${ppn}ppn_${tpr}tpr_${launchID}launchID"
  # 'PreFile' requires NumNodes specification because in the 'Pre' stage, we want to keep the data for different node counts separate.
  local PreFile="${tag1}_${scale}_${matrixDimM}_${matrixDimN}_${curInverseCutOffMult}_${pDimD}_${pDimC}_${ppn}_${tpr}_${NumNodes}nodes"
  local PostFile="${tag1}_${scale}_${matrixDimMorig}_${matrixDimNorig}_${curInverseCutOffMult}_${pDimDorig}_${pDimCorig}_${ppn}_${tpr}"
  local UpdatePlotFile1="${tag1}_${scale}_${matrixDimMorig}_${matrixDimNorig}_${curInverseCutOffMult}_${pDimCorig}"
  local UpdatePlotFile2="${tag1}_${scale}_${matrixDimMorig}_${matrixDimNorig}_${ppn}_${tpr}"

  # Special corner case that only occurs for weak scaling, where invCutOff can increment abruptly to a value its never been before.
  isUniqueTag=1
  collectPlotTagArrayLen=${#collectPlotTags[@]}
  for ((ii=0;ii<${collectPlotTagArrayLen};ii++));
  do
    if [ "${PostFile}" == "${collectPlotTags[${ii}]}" ];
    then
      isUniqueTag=0
    fi
  done
  if [ ${isUniqueTag} -eq 1 ];
  then
    #echo "HERE, plotTags -- ${collectPlotTags[@]}"
    collectPlotTags+=(${PostFile})
  fi

  # Plot instructions only need a single output per scaling study
  if [ ${nodeIndex} == 0 ] || [ ${isUniqueTag} -eq 1 ];
  then
    WriteMethodDataForPlotting 0 ${UpdatePlotFile1} ${UpdatePlotFile2} ${tag1} ${PostFile} ${pDimD} ${pDimC} ${curInverseCutOffMult} ${ppn} ${tpr}
    TemporaryDCplotInfo ${scaleRegime} ${nodeIndex} ${nodeCount} ${pDimDorig} ${pDimCorig} ${WScounterOffset}
    writePlotFileName ${PostFile} $SCRATCH/${testName}/plotInstructions.sh 1  
  fi

  WriteMethodDataForCollectingStage1 ${tag1} ${PreFile} ${PreFile}_perf ${PreFile}_numerics $SCRATCH/${testName}/collectInstructionsStage1.sh
  WriteMethodDataForCollectingStage2 ${launchID} ${tag1} ${PreFile} ${PreFile}_perf ${PreFile}_numerics ${PostFile} ${PostFile}_perf ${PostFile}_numerics $SCRATCH/${testName}/collectInstructionsStage2.sh
  launchJobsPortal ${binaryPath} ${tag1} ${fileString} ${launchID} ${NumNodes} ${ppn} ${tpr} ${matrixDimM} ${matrixDimN} ${bcDim} ${curInverseCutOffMult} 0 ${pDimD} ${pDimC} ${numIterations} $SCRATCH/${testName}/${fileString}
  writePlotFileName ${fileString} $SCRATCH/${testName}/collectInstructionsStage1.sh 0
  curInverseCutOffMult=$(( ${curInverseCutOffMult} + 1 ))
}


# For CFR3D
camfs_cfr3d () {
  # launch CFR3D
  local scale=${1}
  local binaryPath=${2}
  local numIterations=${3}
  local launchID=${4}
  local NumNodes=${5}
  local ppn=${6}
  local tpr=${7}
  local matrixDimM=${8}
  local matrixDimMorig=${9}
  local cubeDimorig=${10}
  local cubeDim=${11}
  local nodeIndex=${12}
  local scaleRegime=${13}
  local nodeCount=${14}
  local bcDim=0
  local tag1="cfr3d"

  curInverseCutOffMult=0
  # Set up the file string that will store the local benchmarking results
  local fileString="DataFiles/results_${tag1}_${scale}_${NumNodes}nodes_${matrixDimM}dimM_${curInverseCutOffMult}inverseCutOffMult_0bcMult_0panelDimMult_${cubeDim}cubeDim_${numIterations}numIter_${ppn}ppn_${tpr}tpr_${curLaunchID}launchID"
  # 'PreFile' requires NumNodes specification because in the 'Pre' stage, we want to keep the data for different node counts separate.
  local PreFile="${tag1}_${scale}_${matrixDimM}_${curInverseCutOffMult}_${cubeDim}_${ppn}_${tpr}_${NumNodes}nodes"
  local PostFile="${tag1}_${scale}_${matrixDimMorig}_${curInverseCutOffMult}_${cubeDimorig}_${ppn}_${tpr}"
  local UpdatePlotFile1="${tag1}_${scale}_${matrixDimMorig}_${curInverseCutOffMult}_${cubeDimorig}"
  local UpdatePlotFile2="${tag1}_${scale}_${matrixDimMorig}_${ppn}_${tpr}"

  # Plot instructions only need a single output per scaling study
  if [ ${nodeIndex} == 0 ];
  then
    WriteMethodDataForPlotting 0 ${UpdatePlotFile1} ${UpdatePlotFile2} ${tag1} ${PostFile} ${cubeDim} ${curInverseCutOffMult} ${ppn} ${tpr}
    writePlotFileName ${PostFile} $SCRATCH/${testName}/plotInstructions.sh 1
  fi

  WriteMethodDataForCollectingStage1 ${tag1} ${PreFile} ${PreFile}_perf ${PreFile}_numerics $SCRATCH/${testName}/collectInstructionsStage1.sh
  WriteMethodDataForCollectingStage2 ${launchID} ${tag1} ${PreFile} ${PreFile}_perf ${PreFile}_numerics ${PostFile} ${PostFile}_perf ${PostFile}_numerics $SCRATCH/${testName}/collectInstructionsStage2.sh
  # Don't pass in 'cubeDim', because this is inferred based on the number of processes, as its just the cube root
  launchJobsPortal ${binaryPath} ${tag1} ${fileString} ${curLaunchID} ${NumNodes} ${ppn} ${tpr} ${matrixDimM} ${bcDim} ${curInverseCutOffMult} 0 ${numIterations} $SCRATCH/${testName}/${fileString}
  writePlotFileName ${fileString} $SCRATCH/${testName}/collectInstructionsStage1.sh 0
  curInverseCutOffMult=$(( ${curInverseCutOffMult} + 1 ))
}


# For MM3D
camfs_mm3d () {
  # launch CFR3D
  local scale=${1}
  local binaryPath=${2}
  local numIterations=${3}
  local launchID=${4}
  local NumNodes=${5}
  local ppn=${6}
  local tpr=${7}
  local gemmORtrmm=${8}
  local algChoice=${9}
  local matrixDimM=${10}
  local matrixDimN=${11}
  local matrixDimK=${12}
  local matrixDimMorig=${13}
  local matrixDimNorig=${14}
  local matrixDimKorig=${15}
  local cubeDimorig=${16}
  local cubeDim=${17}
  local nodeIndex=${18}
  local scaleRegime=${19}
  local nodeCount=${20}
  local tag1="mm3d"

  # Set up the file string that will store the local benchmarking results
  local fileString="DataFiles/results_${tag1}_${scale}_${NumNodes}nodes_${matrixDimM}dimM_${matrixDimN}dimN_${matrixDimK}dimK_${cubeDim}cubeDim_${numIterations}numIter_${ppn}ppn_${tpr}tpr_${curLaunchID}launchID"
  # 'PreFile' requires NumNodes specification because in the 'Pre' stage, we want to keep the data for different node counts separate.
  local PreFile="${tag1}_${scale}_${matrixDimM}_${matrixDimN}_${matrixDimK}_${cubeDim}_${ppn}_${tpr}_${NumNodes}nodes"
  local PostFile="${tag1}_${scale}_${matrixDimMorig}_${matrixDimNorig}_${matrixDimKorig}_${cubeDimorig}_${ppn}_${tpr}"
  local UpdatePlotFile1="${tag1}_${scale}_${matrixDimMorig}_${matrixDimNorig}_${matrixDimKorig}_${cubeDimorig}"
  local UpdatePlotFile2="${tag1}_${scale}_${matrixDimMorig}__${matrixDimNorig}_${matrixDimKorig}${ppn}_${tpr}"

  # Plot instructions only need a single output per scaling study
  if [ ${nodeIndex} == 0 ];
  then
    WriteMethodDataForPlotting 0 ${UpdatePlotFile1} ${UpdatePlotFile2} ${tag1} ${PostFile} ${cubeDim} ${ppn} ${tpr}
    writePlotFileName ${PostFile} $SCRATCH/${testName}/plotInstructions.sh 1
  fi

  WriteMethodDataForCollectingStage1 ${tag1} ${PreFile} ${PreFile}_perf ${PreFile}_numerics $SCRATCH/${testName}/collectInstructionsStage1.sh
  WriteMethodDataForCollectingStage2 ${launchID} ${tag1} ${PreFile} ${PreFile}_perf ${PreFile}_numerics ${PostFile} ${PostFile}_perf ${PostFile}_numerics $SCRATCH/${testName}/collectInstructionsStage2.sh
  # Don't pass in 'cubeDim', because this is inferred based on the number of processes, as its just the cube root
  launchJobsPortal ${binaryPath} ${tag1} ${fileString} ${curLaunchID} ${NumNodes} ${ppn} ${tpr} ${gemmORtrmm} ${algChoice} ${matrixDimM} ${matrixDimN} ${matrixDimK} ${numIterations} $SCRATCH/${testName}/${fileString}
  writePlotFileName ${fileString} $SCRATCH/${testName}/collectInstructionsStage1.sh 0
}


candmc_bsqr () {
  # launch scaLAPACK_QR
  local binaryTag=${1}
  local scale=${2}
  local binaryPath=${3}
  local numIterations=${4}
  local launchID=${5}
  local NumNodes=${6}
  local ppn=${7}
  local tpr=${8}
  local matrixDimM=${9}
  local matrixDimN=${10}
  local matrixDimMorig=${11}
  local matrixDimNorig=${12}
  local numProwsorig=${13}
  local numPcolsorig=${14}
  local numProws=${15}
  local minBlockSize=${16}
  local maxBlockSize=${17}
  local nodeIndex=${18}
  local scaleRegime=${19}
  local nodeCount=${20}
  local tag1="bsqr"

  for ((k=${minBlockSize}; k<=${maxBlockSize}; k*=2))
  do
    # Set up the file string that will store the local benchmarking results
    local fileString="DataFiles/results_${tag1}_${scale}_${NumNodes}nodes_${matrixDimM}dimM_${matrixDimN}dimN_${numProws}numProws_${k}bSize_${numIterations}numIter_${ppn}ppn_${tpr}tpr_${curLaunchID}launchID"
    # 'PreFile' requires NumNodes specification because in the 'Pre' stage, we want to keep the data for different node counts separate.
    local PreFile="${tag1}_${scale}_${matrixDimM}_${matrixDimN}_${numProws}_${k}_${ppn}_${tpr}_${NumNodes}nodes"
    local PostFile="${tag1}_${scale}_${matrixDimMorig}_${matrixDimNorig}_${numProwsorig}_${k}_${ppn}_${tpr}"
    local UpdatePlotFile1="${tag1}_${scale}_${matrixDimMorig}_${matrixDimNorig}_${numPcolsorig}_${k}"
    local UpdatePlotFile2="${tag1}_${scale}_${matrixDimMorig}_${matrixDimNorig}_${ppn}_${tpr}"

    # Plot instructions only need a single output per scaling study
    if [ ${nodeIndex} == 0 ];
    then
      WriteMethodDataForPlotting 0 ${UpdatePlotFile1} ${UpdatePlotFile2} ${tag1} ${PostFile} ${numProws} ${k} ${ppn} ${tpr}
      writePlotFileNameScalapackQR ${PostFile} $SCRATCH/${testName}/plotInstructions.sh 1
    fi

    WriteMethodDataForCollectingStage1 ${tag1} ${PreFile} ${PreFile}_NoFormQ ${PreFile}_FormQ $SCRATCH/${testName}/collectInstructionsStage1.sh
    WriteMethodDataForCollectingStage2 ${launchID} ${tag1} ${PreFile} ${PreFile}_NoFormQ ${PreFile}_FormQ ${PostFile} ${PostFile}_NoFormQ ${PostFile}_FormQ $SCRATCH/${testName}/collectInstructionsStage2.sh
    launchJobsPortal ${binaryPath} ${tag1} ${fileString} ${curLaunchID} ${NumNodes} ${ppn} ${tpr} ${matrixDimM} ${matrixDimN} ${k} ${numIterations} 0 ${numProws} 1 0 $SCRATCH/${testName}/${fileString}
    writePlotFileNameScalapackQR ${fileString} $SCRATCH/${testName}/collectInstructionsStage1.sh 0
  done
}

candmc_bscf () {
  # launch scaLAPACK_CF
  local scale=${1}
  local binaryPath=${2}
  local numIterations=${3}
  local launchID=${4}
  local NumNodes=${5}
  local ppn=${6}
  local tpr=${7}
  local matrixDimM=${8}
  local matrixDimMorig=${8}
  local minBlockSize=${9}
  local maxBlockSize=${10}
  local nodeIndex=${11}
  local scaleRegime=${12}
  local nodeCount=${13}
  local tag1="bscf"

  for ((k=${minBlockSize}; k<=${maxBlockSize}; k*=2))
  do
    # Set up the file string that will store the local benchmarking results
    local fileString="DataFiles/results_${tag1}_${scale}_${NumNodes}nodes_${matrixDimM}dimM_${k}bSize_${numIterations}numIter_${ppn}ppn_${tpr}tpr_${curLaunchID}launchID"
    # 'PreFile' requires NumNodes specification because in the 'Pre' stage, we want to keep the data for different node counts separate.
    local PreFile="${tag1}_${scale}_${matrixDimM}_${k}_${ppn}_${tpr}_${NumNodes}nodes"
    local PostFile="${tag1}_${scale}_${matrixDimMorig}_${k}_${ppn}_${tpr}"
    local UpdatePlotFile1="${tag1}_${scale}_${matrixDimMorig}_${k}"
    local UpdatePlotFile2="${tag1}_${scale}_${ppn}_${tpr}"

    if [ ${nodeIndex} == 0 ];
    then
      # Write to plotInstructions file
      WriteMethodDataForPlotting 0 ${UpdatePlotFile1} ${UpdatePlotFile2} ${tag1} ${PostFile} ${k} ${ppn} ${tpr}
      writePlotFileNameScalapackCholesky ${PostFile} $SCRATCH/${testName}/plotInstructions.sh 1
    fi

    WriteMethodDataForCollectingStage1 ${tag1} ${PreFile} ${PreFile} ${PreFile}_blah $SCRATCH/${testName}/collectInstructionsStage1.sh
    WriteMethodDataForCollectingStage2 ${launchID} ${tag1} ${PreFile} ${PreFile} ${PreFile}_blah ${PostFile} ${PostFile} ${PostFile} $SCRATCH/${testName}/collectInstructionsStage2.sh
    launchJobsPortal ${binaryPath} ${tag1} ${fileString} ${curLaunchID} ${NumNodes} ${ppn} ${tpr} ${matrixDimM} ${k} ${numIterations} $SCRATCH/${testName}/${fileString}
    writePlotFileNameScalapackCholesky ${fileString} $SCRATCH/${testName}/collectInstructionsStage1.sh 0
  done
}
