# For CFR3D
launchcfr3d () {
  # launch CFR3D
  local scale=\${1}
  local binaryPath=\${2}
  local numIterations=\${3}
  local launchID=\${4}
  local NumNodes=\${5}
  local ppn=\${6}
  local tpr=\${7}
  local matrixDimM=\${8}
  local matrixDimMorig=\${9}
  local cubeDimorig=\${10}
  local cubeDim=\${11}
  local nodeIndex=\${12}
  local scaleRegime=\${13}
  local nodeCount=\${14}
  local bcDim=0

  # Next: Based on pDimC, decide on invCutOff parameter, which will range from 0 to a max of 2 for now
  invCutOffLoopMax=0
  if [ \${cubeDim} -le 2 ];
  then
    invCutOffLoopMax=0
  elif [ \${cubeDim} -eq 4 ];
  then
    invCutOffLoopMax=1
  else
    invCutOffLoopMax=2
    #invCutOffLoopMax=\$(( \${cubeDim} / 2 ))
    #invCutOffLoopMax=\$(( \${invCutOffLoopMax} - 1 ))
  fi

  curInverseCutOffMult=0
  while [ \${curInverseCutOffMult} -le \${invCutOffLoopMax} ];
  do
    # Set up the file string that will store the local benchmarking results
    local fileString="DataFiles/results_${tag3}_\${scale}_\${NumNodes}nodes_\${matrixDimM}dimM_\${curInverseCutOffMult}inverseCutOffMult_0bcMult_0panelDimMult_\${cubeDim}cubeDim_\${numIterations}numIter_\${ppn}ppn_\${tpr}tpr_\${curLaunchID}launchID"
    # 'PreFile' requires NumNodes specification because in the 'Pre' stage, we want to keep the data for different node counts separate.
    local PreFile="${tag3}_\${scale}_\${matrixDimM}_\${curInverseCutOffMult}_\${cubeDim}_\${ppn}_\${tpr}_\${NumNodes}nodes"
    local PostFile="${tag3}_\${scale}_\${matrixDimMorig}_\${curInverseCutOffMult}_\${cubeDimorig}_\${ppn}_\${tpr}"
    local UpdatePlotFile1="${tag3}_\${scale}_\${matrixDimMorig}_\${curInverseCutOffMult}_\${cubeDimorig}"
    local UpdatePlotFile2="${tag3}_\${scale}_\${matrixDimMorig}_\${ppn}_\${tpr}"

    # Plot instructions only need a single output per scaling study
    if [ \${nodeIndex} == 0 ];
    then
      WriteMethodDataForPlotting 0 \${UpdatePlotFile1} \${UpdatePlotFile2} ${tag3} \${PostFile} \${cubeDim} \${curInverseCutOffMult} \${ppn} \${tpr}
      writePlotFileName \${PostFile} $SCRATCH/${fileName}/plotInstructions.sh 1
    fi

    WriteMethodDataForCollectingStage1 ${tag3} \${PreFile} \${PreFile}_perf \${PreFile}_numerics $SCRATCH/${fileName}/collectInstructionsStage1.sh
    WriteMethodDataForCollectingStage2 \${launchID} ${tag3} \${PreFile} \${PreFile}_perf \${PreFile}_numerics \${PostFile} \${PostFile}_perf \${PostFile}_numerics $SCRATCH/${fileName}/collectInstructionsStage2.sh
    # Don't pass in 'cubeDim', because this is inferred based on the number of processes, as its just the cube root
    launchJobsPortal \${binaryPath} ${tag3} \${fileString} \${curLaunchID} \${NumNodes} \${ppn} \${tpr} \${matrixDimM} \${bcDim} \${curInverseCutOffMult} 0 \${numIterations} $SCRATCH/${fileName}/\${fileString}
    writePlotFileName \${fileString} $SCRATCH/${fileName}/collectInstructionsStage1.sh 0
    curInverseCutOffMult=\$(( \${curInverseCutOffMult} + 1 ))
  done
}
