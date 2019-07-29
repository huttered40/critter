import os
from subprocess import call


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
  invCutOffLoopMax=0
  # Note: I am just not seeing enough performance boost at the moment to warrant launching the varying invCutOff runs, especially for Critter
  if [ "${tuneInvCutOff}" == "y" ];
  then
    if [ ${pDimC} -le 2 ];
    then
      invCutOffLoopMax=0
    elif [ ${pDimC} -eq 4 ];
    then
      invCutOffLoopMax=1
    else
      invCutOffLoopMax=2
      #invCutOffLoopMax=$(( ${pDimC} / 2 ))
      #invCutOffLoopMax=$(( ${invCutOffLoopMax} - 1 ))
    fi
  fi

  curInverseCutOffMult=0
  while [ ${curInverseCutOffMult} -le ${invCutOffLoopMax} ];
  do
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
  done
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

  # Next: Based on pDimC, decide on invCutOff parameter, which will range from 0 to a max of 2 for now
  invCutOffLoopMax=0
  if [ ${cubeDim} -le 2 ];
  then
    invCutOffLoopMax=0
  elif [ ${cubeDim} -eq 4 ];
  then
    invCutOffLoopMax=1
  else
    invCutOffLoopMax=2
    #invCutOffLoopMax=$(( ${cubeDim} / 2 ))
    #invCutOffLoopMax=$(( ${invCutOffLoopMax} - 1 ))
  fi

  curInverseCutOffMult=0
  while [ ${curInverseCutOffMult} -le ${invCutOffLoopMax} ];
  do
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
  done
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



class bench(object):
    """
    """
    def __init__(self,\
                 CritterPath,\
		 MachineType,\
		 LibraryTypeList,\
		 fileID,\
		 roundID,\
		 NumLaunchesPerBinary,\
		 numTests,\
		 numHours,\
		 numMinutes,\
		 numSeconds,\
		 email,\
		 dataType,\
		 intType,\
		 analyzeDecision1,\
		 analyzeDecision2,\
		 mpiType,\
		 minPEcountPerNode,\
		 maxPEcountPerNode,\
		 nodeMinList,\
		 nodeMaxList,\
		 ppnMinList,\
		 ppnMaxList,\
		 tprMinList,\
		 tprMaxList,\
		 nodeScaleFactorList,\
		 ppnScaleFactorList,\
		 tprScaleFactorList,\
		 SubmitToQueue,\
		 AlgorithmList):
       """
       CritterPath - specify full path to Critter repository
                   - do not include a '/' after 'critter'
                   - specified as a string

       MachineType - machine type specified in ../experiments/machines/

       LibraryTypeList - list of types specified in ../experiments/libraries/

       fileID - base name of the directory inside which all data/scripts will be stored
              - will appear inside the SCRATCH directory
              - specified as a string

       roundID - set to '1' unless performing piecewise testing (launching same job separately) to enhance performance reproducibility

       NumLaunchesPerBinary - set to '1' unless performing testing to enhance performance reproducibility (by launching same jobs multiple times)
                            - different from 'roundID' because the former did not launch at the same time, but waited via a separate launch of bench.sh

       numTests - number of scaling studies
                - for example, weak scaling and strong scaling, even across the same variants, constitute separate tests

       numHours,numMinutes,numSeconds - time for job
                                      - must be specified as two-digit strings

       email - specify email address that you'd like job updates to appear

       dataType - float[0], double[1], complex<float>[2], complex<double>[3]
                - only relevant if test file uses an environment variable to specify this

       intType - int[0], int64_t[1]
               - only relevant if test file uses an environment variable to specify this

       analyzeDecision1 - profile using Critter

       analyzeDecision2 - profile using TAU
                        - not currently supported

       mpiType - specify 'mpi' (unless on Porter, then 'ampi' is available)

       minPEcountPerNode/maxPEcountPerNode - specify the min and max number of processing elements (processes per node x threads per process)

       nodeMinList,nodeMaxList - min/max number of nodes for each test

       ppnMinList,ppnMaxList - min/max number of processes-per-node for each node count for each test

       tprMinList,tprMaxList - min/max number of threads-per-process for each node count for each test

       nodeScaleFactorList - scaling factor to apply to the number of nodes for each test

       ppnScaleFactorList - scaling factor to apply to the number of MPI processes per node (ppn) for each test

       tprScaleFactorList - scaling factor to apply to the number of threads per MPI rank (tpr) for each test

       SubmitToQueue - '1' to submit jobs to queue, '0' to not submit to queue

       AlgorithmList - list of lists of lists of inputs, must be of length 'numTests'
                          - each inner list holds the inputs for the corresponding algorithm in AlgorithmList
       """
       self.CritterPath = CritterPath
       self.MachineType = MachineType
       self.LibraryTypeList = LibraryTypeList
       self.fileID = fileID
       self.roundID = roundID
       self.NumLaunchesPerBinary = NumLaunchesPerBinary
       self.numTests = numTests
       self.numHours = numHours
       self.numMinutes = numMinutes
       self.numSeconds = numSeconds
       self.email = email
       self.dataType = dataType
       self.intType = intType
       self.analyzeDecision1 = analyzeDecision1
       self.analyzeDecision2 = analyzeDecision2
       self.mpiType = mpiType
       self.minPEcountPerNode = minPEcountPerNode
       self.maxPEcountPerNode = maxPEcountPerNode
       self.nodeMinList = nodeMinList
       self.nodeMaxList = nodeMaxList
       self.ppnMinList = ppnMinList
       self.ppnMaxList = ppnMaxList
       self.tprMinList = tprMinList
       self.tprMaxList = tprMaxList
       self.nodeScaleFactorList = nodeScaleFactorList
       self.ppnScaleFactorList = ppnScaleFactorList
       self.tprScaleFactorList = tprScaleFactorList
       self.SubmitToQueue = SubmitToQueue
       self.AlgorithmList = AlgorithmList
       dateStr=$(date +%Y-%m-%d-%H_%M_%S)
       self.testName="%s_%s_%s_round%d"%(fileID,dateStr,self.MachineType.machineName,roundID)
       self.testNameAllRounds="%s_%s"%(fileID,self.MachineType.machineName)

    def __portal(self,func,op):
        """
        """
        # Submit all scripts
        .. due to new format, need to watch out for repeats
        curLaunchID=1
        while [ ${curLaunchID} -le ${NumLaunchesPerBinary} ];
            curNumNodes=${minNumNodes}
            listIndex=0
            while [ ${curNumNodes} -le ${maxNumNodes} ];
                curPPN=${ppnMinList[${listIndex}]}
                ppnMax=${ppnMaxList[${listIndex}]}
                while [ ${curPPN} -le ${ppnMax} ];
                    curTPR=${tprMinList[${listIndex}]}
                    tprMax=${tprMaxList[${listIndex}]}
                    while [ ${curTPR} -le ${tprMax} ];
                        # Make sure we are in a suitable range
                        numPEsPerNode=$(( ${curPPN} * ${curTPR} ))
                        if [ ${minPEcountPerNode} -le ${numPEsPerNode} ] && [ ${maxPEcountPerNode} -ge ${numPEsPerNode} ];

                            .. user 'op' to differentiate between whether to pass 9 args (to op=2) or not
                            .. note that we will want to append to the script files, and this is important, since python might have a special tag for that

                            FullScriptName=${testName}/script_${fileID}id_${roundID}round_${curLaunchID}launchID_${curNumNodes}nodes_${curPPN}ppn_${curTPR}tpr.${BatchFileExtension}
                            scriptName=$SCRATCH/${testName}/script_${fileID}id_${roundID}round_\${curLaunchID}launchID_\${curNumNodes}nodes_\${curPPN}ppn_\${curTPR}tpr.${BatchFileExtension}
                            def script(scriptFile,testName,curNumNodes,curPPN,curTPR,numPEsPerNode,numHours,numMinutes,numSeconds):
                        curTPR=$(( ${curTPR} * ${tprScaleFactor} ))
                    curPPN=$(( ${curPPN} * ${ppnScaleFactor} ))
                curNumNodes=$(( ${curNumNodes} * ${nodeScaleFactor} ))
                listIndex=$(( ${listIndex} + 1 ))
            curLaunchID=$(( ${curLaunchID} + 1 ))

    def queue_submit(self):
        """
        """
	if (self.SubmitToQueue == 1):
	  call("mkdir %s/%s/bin"%(os.environ["SCRATCH"],self.testName))
	  call("mv ../Tests/%s/* %s/%s/bin"%(self.testName,os.environ["SCRATCH"],self.testName))
          portal(self.MachineType.queue)

    def launch(self):
        """
        """
        self.MachineType.set()

        if (self.mpiType == "mpi"):
            os.environ["MPITYPE"] = "MPI_TYPE"
        elif (self.mpiType == "ampi"):
            os.environ["MPITYPE"] = "AMPI_TYPE"

        call("mkdir ../Tests/%s"%(self.testName))
        call("mkdir ../Tests/%s/bin"%(self.testName))

        if (dataType == 0):
            os.environ["DATATYPE"] = "FLOAT_TYPE"
        elif (dataType == 1)
            os.environ["DATATYPE"] = "DOUBLE_TYPE"
        elif (dataType == 2):
            os.environ["DATATYPE"] = "COMPLEX_FLOAT_TYPE"
        elif (dataType == 3):
            os.environ["DATATYPE"] = "COMPLEX_DOUBLE_TYPE"
        if (intType == 0):
            os.environ["INTTYPE"] = "INT_TYPE"
        elif (intType == 1):
            os.environ["INTTYPE"] = "INT64_T_TYPE"

        for lib in self.LibraryTypeList:
            os.environ["PROFTYPE"]="PERFORMANCE"
	    profType="P"
            # export SPECIAL_SCALA_ARG=REF
            lib.build()
            if [ self.analyzeDecision1 == 1 ];
            then
                profType="PC"
                os.environ["PROFTYPE"]="CRITTER"
                lib.build()
            fi

            if [ ${analyzeDecision2} == 1 ];
            then
                profType=profType+"T"
                export PROFTYPE=PROFILE
                os.environ["PROFTYPE"]="PROFILE"
                lib.build()
            fi

        os.environ["BINARYPATH"] = os.environ["SCRATCH"] + "/%s/bin/"%(self.testName)

        # collectData.sh will always be a single line, just a necessary intermediate step.
        ## echo "bash $SCRATCH/${testName}/collectInstructionsStage1.sh | bash PackageDataRemoteStage1.sh" > collectData.sh

        # Launch the generated script
        #bash $SCRATCH/${testName}.sh

        call("mkdir %s/%s/"%(os.environ["SCRATCH"],self.testName),shell=True)
        call("mkdir %s/%s/DataFiles/"%(os.environ["SCRATCH"],self.testName),shell=True)
        portal(self.MachineType.script)

        .. need to open 3 files for appending: plotInstructions.sh, collectInstructionsStage1.sh, collectInstructionsStage2.sh

        # Note: in future, I may want to decouple numBinaries and numPlotTargets, but only when I find it necessary
        # Write to Plot Instructions file, for use by SCAPLOT makefile generator
        echo "echo \"1\"" > $SCRATCH/${testName}/plotInstructions.sh
        echo "echo \"${testNameAllRounds}\"" >> $SCRATCH/${testName}/plotInstructions.sh
        echo "echo \"${numTests}\"" >> $SCRATCH/${testName}/plotInstructions.sh
        echo "echo \"${machineName}\"" >> $SCRATCH/${testName}/plotInstructions.sh
        echo "echo \"${profType}\"" >> $SCRATCH/${testName}/plotInstructions.sh
        echo "echo \"${nodeScaleFactor}\"" >> $SCRATCH/${testName}/plotInstructions.sh

        WriteHeaderForCollection $SCRATCH/${testName}/collectInstructionsStage1.sh
        WriteHeaderForCollection $SCRATCH/${testName}/collectInstructionsStage2.sh

        for test in range(1,numTests+1):
            print("Test %d\n"%(i))

            # Nodes
            startNumNodes,endNumNodes comes from the corresponding test

            scaleRegime,scale comes from algList .. problem here is that the .sh files want them as test-level, not alg-level

            echo "echo \"\${scale}\"" >> $SCRATCH/${testName}/plotInstructions.sh

            nodeCount= get from len of list
            echo "echo \"\${nodeCount}\" " >> $SCRATCH/${testName}/plotInstructions.sh

            curNumNodes=\${startNumNodes}
            for ((j=0; j<\${nodeCount}; j++))
                echo "echo \"\${curNumNodes}\"" >> $SCRATCH/${testName}/plotInstructions.sh
                curNumNodes=\$(( \${curNumNodes} * ${nodeScaleFactor} ))

            matrixDimM,matrixDimN,numIterations, comes from AlgList

            echo "echo \"\${matrixDimM}\"" >> $SCRATCH/${testName}/plotInstructions.sh
            echo "echo \"\${matrixDimN}\"" >> $SCRATCH/${testName}/plotInstructions.sh


            .. basically the loop below should go over all algorithms, BUT it should also go over all starting input parameters, which is set by the Alg
            .. need a separate loop for that
            .. also need a check for the last input parameters. maybe just a while (), then grab
            j=1
            while [ 1 -eq 1 ];		# Loop iterates until user says stop
                echo -e "\nStage #\${j}"

                # Echo for SCAPLOT makefile generator
                binaryTag=self.AlgorithmList[...].Tag
                echo "echo \"\${binaryTag}\"" >> $SCRATCH/${testName}/collectInstructionsStage1.sh
                echo "echo \"\${binaryTag}\"" >> $SCRATCH/${testName}/collectInstructionsStage2.sh
                echo "echo \"\${binaryTag}\"" >> $SCRATCH/${testName}/plotInstructions.sh

                ..binaryPath=${BINARYPATH}\${binaryTag}_${machineName}
                if [ "${machineName}" == "PORTER" ];
                    binaryPath=\${binaryPath}_${mpiType}
                elif [ "${machineName}" == "BLUEWATERS" ];
                    # special case, only for CAMFS, not for bench_scalapack routines
                    if [ "\${binaryTag}" == "camfs_cacqr2" ] || [ "\${binaryTag}" == "camfs_cfr3d" ];
                        binaryPath=${BINARYPATH}\${binaryTag}_${machineName}_${GPU}

                .. can we re-use portal here for the outer-loop-structure?
                for ((curLaunchID=1; curLaunchID<=${NumLaunchesPerBinary}; curLaunchID+=1));
                    # Initialize all possible variables that change with node count
                    # shared
                    nodeIndex=0
                    curMatrixDimM=\${matrixDimM}
                    curMatrixDimN=\${matrixDimN}
                    if [ \${binaryTag} == 'mm3d' ];
                        curMatrixDimK=\${matrixDimK}
                    # cacqr2
                    pDimCArray=()
                    pDimCArrayOrig=()
                    rangePdimClen=0
                    if [ \${binaryTag} == 'camfs_cacqr2' ];
                        for ((w=\${startStartPdimC}; w<=\${endStartPdimC}; w*=2));
                            pDimCArray+=(\${w})
                            pDimCArrayOrig+=(\${w})
                            rangePdimClen=\$(( \${rangePdimClen} + 1 ))
                    # bsqr/rsqr
                    numPcolsArray=()
                    numPcolsArrayOrig=()
                    rangeNumPcolslen=0
                    if [ \${binaryTag} == 'bsqr' ] || [ \${binaryTag} == 'rsqr' ];
                        for ((w=\${startStartNumPcols}; w<=\${endStartNumPcols}; w*=2));
                            numPcolsArray+=(\${w})
                            numPcolsArrayOrig+=(\${w})
                            rangeNumPcolslen=\$(( \${rangeNumPcolslen} + 1 ))
                    # cfr3d
                    curCubeDim=\${cubeDim}
                    for ((curNumNodes=\${startNumNodes}; curNumNodes<=\${endNumNodes}; curNumNodes*=${nodeScaleFactor}));
                        minPPN=\${ppnMinListRunTime[\${nodeIndex}]}
                        maxPPN=\${ppnMaxListRunTime[\${nodeIndex}]}
                        for ((curPPN=\${minPPN}; curPPN<=\${maxPPN}; curPPN*=${ppnScaleFactor}));
                            numProcesses=\$(( \${curNumNodes} * \${curPPN} ))
                            StartingNumProcesses=\$(( \${startNumNodes} * \${curPPN} ))

                            minTPR=\${tprMinListRunTime[\${nodeIndex}]}
                            maxTPR=\${tprMaxListRunTime[\${nodeIndex}]}
                            for ((curTPR=\${minTPR}; curTPR<=\${maxTPR}; curTPR*=${tprScaleFactor}));
                                # Make sure we are in a suitable range
                                numPEsPerNode=\$(( \${curPPN} * \${curTPR} ))
                                if [ ${minPEcountPerNode} -le \${numPEsPerNode} ] && [ ${maxPEcountPerNode} -ge \${numPEsPerNode} ];
                                    if [ \${binaryTag} == 'camfs_cacqr2' ];
                                        # Below: note that the STARTING dimC is being changed. The parameters that aren't solely dependent on the node count are
                                        #   changed here and not in launchTag***
                                        for ((w=0; w<\${rangePdimClen}; w+=1));
                                            pDimC=\${pDimCArray[\${w}]}
                                            pDimCsquared=\$(( \${pDimC} * \${pDimC} ))
                                            pDimD=\$(( \${numProcesses} / \${pDimCsquared} ))

                                            # Special check because performance for 16 PPN, 4 TPR shows superior performance for the skinniest grid
                                            isSpecial=1

                                            # Check if pDimC is too big. If so, pDimD will be 0
                                            if [ \${pDimD} -ge \${pDimC} ] && [ \${isSpecial} == 1 ];
                                                originalPdimC=\${pDimCArrayOrig[\${w}]}
                                                originalPdimCsquared=\$(( \${originalPdimC} * \${originalPdimC} ))
                                                originalPdDimD=\$(( \${StartingNumProcesses} / \${originalPdimCsquared} ))
                                                \${binaryTag} \${scale..} \${binaryPath} \${numIterations} \${curLaunchID} \${curNumNodes} \${curPPN} \${curTPR} \${curMatrixDimM} \${curMatrixDimN} \${matrixDimM} \${matrixDimN} \${originalPdDimD} \${originalPdimC} \${pDimD} \${pDimC} \${nodeIndex} \${scaleRegime..} \${nodeCount} \${WShelpcounter} \${invCutOffDec}
                                    elif [ \${binaryTag} == 'bsqr' ] || [ \${binaryTag} == 'rsqr' ];
                                        # Special case to watch out for.
                                        for ((w=0; w<\${rangeNumPcolslen}; w+=1));
                                        numPcols=\${numPcolsArray[\${w}]}
                                        numProws=\$(( \${numProcesses} / \${numPcols} ))
                                        isSpecial=1
                                        if [ \${numPcols} -le \${numProws} ] && [ \${isSpecial} == 1 ];
                                            originalNumPcols=\${numPcolsArrayOrig[\${w}]}
                                            originalNumProws=\$(( \${StartingNumProcesses} / \${originalNumPcols} ))
                                            sharedBinaryTag="candmc_bsqr"	# Even if rsqr, use bsqr and then have the corresponding method use the new argument for binaryTag
                                            \${sharedBinaryTag} \${binaryTag} \${scale..} \${binaryPath} \${numIterations} \${curLaunchID} \${curNumNodes} \${curPPN} \${curTPR} \${curMatrixDimM} \${curMatrixDimN} \${matrixDimM} \${matrixDimN} \${originalNumProws} \${originalNumPcols} \${numProws} \${minBlockSize} \${maxBlockSize} \${nodeIndex} \${scaleRegime..} \${nodeCount}
                                    elif [ \${binaryTag} == 'cfr3d' ];
                                        \${binaryTag} \${scale..} \${binaryPath} \${numIterations} \${curLaunchID} \${curNumNodes} \${curPPN} \${curTPR} \${curMatrixDimM} \${matrixDimM} \${cubeDim} \${curCubeDim} \${nodeIndex} \${scaleRegime..} \${nodeCount}
                                    elif [ \${binaryTag} == 'bscf' ];
                                        \${binaryTag} \${scale..} \${binaryPath} \${numIterations} \${curLaunchID} \${curNumNodes} \${curPPN} \${curTPR} \${curMatrixDimM} \${matrixDimM} \${minBlockSize} \${maxBlockSize} \${nodeIndex} \${scaleRegime..} \${nodeCount}
                                    elif [ \${binaryTag} == 'mm3d' ];
                                        \${binaryTag} \${scale..} \${binaryPath} \${numIterations} \${curLaunchID} \${curNumNodes} \${curPPN} \${curTPR} \${gemmORtrmmChoice} \${bcastORallgatherChoice} \${curMatrixDimM} \${curMatrixDimN} \${curMatrixDimK} \${matrixDimM} \${matrixDimN} \${matrixDimK} \${cubeDim} \${curCubeDim} \${nodeIndex} \${scaleRegime..} \${nodeCount}
                        alg.scale(alg.IndirectIndexFunc(scaleCount))
	                nodeIndex=\$(( \${nodeIndex} + 1 ))
            j=\$(( \${j} + 1 ))
            echo "echo \"1\"" >> $SCRATCH/${testName}/collectInstructionsStage1.sh	# Signals end of the data files for this specific methodID
            echo "echo \"1\"" >> $SCRATCH/${testName}/collectInstructionsStage2.sh	# Signals end of the data files for this specific methodID
            echo "echo \"1\"" >> $SCRATCH/${testName}/plotInstructions.sh	# Signals end of the data files for this specific methodID

        queue_submit()

