# For MM3D
mm3d () {
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
