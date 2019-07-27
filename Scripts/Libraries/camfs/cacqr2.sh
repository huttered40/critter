# For CA-CQR2
launchcqr2 () {
  # launch CQR2
  local scale=\${1}
  local binaryPath=\${2}
  local numIterations=\${3}
  local launchID=\${4}
  local NumNodes=\${5}
  local ppn=\${6}
  local tpr=\${7}
  local matrixDimM=\${8}
  local matrixDimN=\${9}
  local matrixDimMorig=\${10}
  local matrixDimNorig=\${11}
  local pDimDorig=\${12}
  local pDimCorig=\${13}
  local pDimD=\${14}
  local pDimC=\${15}
  local nodeIndex=\${16}
  local scaleRegime=\${17}
  local nodeCount=\${18}
  local WScounterOffset=\${19}
  local tuneInvCutOff=\${20}
  local bcDim=0

  # Next: Based on pDimC, decide on invCutOff parameter, which will range from 0 to a max of 2 for now
  invCutOffLoopMax=0
  # Note: I am just not seeing enough performance boost at the moment to warrant launching the varying invCutOff runs, especially for Critter
  if [ "\${tuneInvCutOff}" == "y" ];
  then
    if [ \${pDimC} -le 2 ];
    then
      invCutOffLoopMax=0
    elif [ \${pDimC} -eq 4 ];
    then
      invCutOffLoopMax=1
    else
      invCutOffLoopMax=2
      #invCutOffLoopMax=\$(( \${pDimC} / 2 ))
      #invCutOffLoopMax=\$(( \${invCutOffLoopMax} - 1 ))
    fi
  fi

  curInverseCutOffMult=0
  while [ \${curInverseCutOffMult} -le \${invCutOffLoopMax} ];
  do
    # Set up the file string that will store the local benchmarking results
    local fileString="DataFiles/results_${tag1}_\${scale}_\${NumNodes}nodes_\${matrixDimM}dimM_\${matrixDimN}dimN_\${curInverseCutOffMult}inverseCutOffMult_0bcMult_0panelDimMult_\${pDimD}pDimD_\${pDimC}pDimC_\${numIterations}numIter_\${ppn}ppn_\${tpr}tpr_\${launchID}launchID"
    # 'PreFile' requires NumNodes specification because in the 'Pre' stage, we want to keep the data for different node counts separate.
    local PreFile="${tag1}_\${scale}_\${matrixDimM}_\${matrixDimN}_\${curInverseCutOffMult}_\${pDimD}_\${pDimC}_\${ppn}_\${tpr}_\${NumNodes}nodes"
    local PostFile="${tag1}_\${scale}_\${matrixDimMorig}_\${matrixDimNorig}_\${curInverseCutOffMult}_\${pDimDorig}_\${pDimCorig}_\${ppn}_\${tpr}"
    local UpdatePlotFile1="${tag1}_\${scale}_\${matrixDimMorig}_\${matrixDimNorig}_\${curInverseCutOffMult}_\${pDimCorig}"
    local UpdatePlotFile2="${tag1}_\${scale}_\${matrixDimMorig}_\${matrixDimNorig}_\${ppn}_\${tpr}"

    # Special corner case that only occurs for weak scaling, where invCutOff can increment abruptly to a value its never been before.
    isUniqueTag=1
    collectPlotTagArrayLen=\${#collectPlotTags[@]}
    for ((ii=0;ii<\${collectPlotTagArrayLen};ii++));
    do
      if [ "\${PostFile}" == "\${collectPlotTags[\${ii}]}" ];
      then
        isUniqueTag=0
      fi
    done
    if [ \${isUniqueTag} -eq 1 ];
    then
      #echo "HERE, plotTags -- \${collectPlotTags[@]}"
      collectPlotTags+=(\${PostFile})
    fi

    # Plot instructions only need a single output per scaling study
    if [ \${nodeIndex} == 0 ] || [ \${isUniqueTag} -eq 1 ];
    then
      WriteMethodDataForPlotting 0 \${UpdatePlotFile1} \${UpdatePlotFile2} ${tag1} \${PostFile} \${pDimD} \${pDimC} \${curInverseCutOffMult} \${ppn} \${tpr}
      TemporaryDCplotInfo \${scaleRegime} \${nodeIndex} \${nodeCount} \${pDimDorig} \${pDimCorig} \${WScounterOffset}
      writePlotFileName \${PostFile} $SCRATCH/${fileName}/plotInstructions.sh 1  
    fi

    WriteMethodDataForCollectingStage1 ${tag1} \${PreFile} \${PreFile}_perf \${PreFile}_numerics $SCRATCH/${fileName}/collectInstructionsStage1.sh
    WriteMethodDataForCollectingStage2 \${launchID} ${tag1} \${PreFile} \${PreFile}_perf \${PreFile}_numerics \${PostFile} \${PostFile}_perf \${PostFile}_numerics $SCRATCH/${fileName}/collectInstructionsStage2.sh
    launchJobsPortal \${binaryPath} ${tag1} \${fileString} \${launchID} \${NumNodes} \${ppn} \${tpr} \${matrixDimM} \${matrixDimN} \${bcDim} \${curInverseCutOffMult} 0 \${pDimD} \${pDimC} \${numIterations} $SCRATCH/${fileName}/\${fileString}
    writePlotFileName \${fileString} $SCRATCH/${fileName}/collectInstructionsStage1.sh 0
    curInverseCutOffMult=\$(( \${curInverseCutOffMult} + 1 ))
  done
}
