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
                 nodeScaleOperatorList,\
                 ppnScaleOperatorList,\
                 tprScaleOperatorList,\
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

        nodeScaleOperatorList - scaling operator to apply to the number of nodes for each test

        ppnScaleOperatorList - scaling operator to apply to the number of MPI processes per node (ppn) for each test

        tprScaleOperatorList - scaling operator to apply to the number of threads per MPI rank (tpr) for each test

        SubmitToQueue - '1' to submit jobs to queue, '0' to not submit to queue

        AlgorithmList - list of lists of algorithm class instances
                      - outer list must be of length 'numTests'
                      - each inner list holds algorithm class instances in a list, and a string specifying a tag as to what kind of scaling is occuring
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
        self.nodeScaleOperatorList=nodeScaleOperatorList
        self.ppnScaleOperatorList=ppnScaleOperatorList
        self.tprScaleOperatorList=tprScaleOperatorList
        self.SubmitToQueue = SubmitToQueue
        self.AlgorithmList = AlgorithmList
        dateStr=$(date +%Y-%m-%d-%H_%M_%S)
        self.testName="%s_%s_%s_round%d"%(fileID,dateStr,self.MachineType.machineName,roundID)
        self.testNameAllRounds="%s_%s"%(fileID,self.MachineType.machineName)

        if (self.mpiType == "mpi"):
            os.environ["MPITYPE"] = "MPI_TYPE"
        elif (self.mpiType == "ampi"):
            os.environ["MPITYPE"] = "AMPI_TYPE"

        # I think these directories serve mainly as a intermediate place to put the binaries
	#   before being moved to SCRATCH
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

        self.MachineType.set()
        os.environ["BINARYPATH"] = os.environ["SCRATCH"] + "/%s/bin/"%(self.testName)

        call("mkdir %s/%s/"%(os.environ["SCRATCH"],self.testName),shell=True)
        call("mkdir %s/%s/DataFiles/"%(os.environ["SCRATCH"],self.testName),shell=True)

        self.PlotInstructionsFile = open("%s/%s/plotInstructions.sh"%(os.environ["SCRATCH"],self.testName),"a+")
        self.CollectInstructions1File = open("%s/%s/plotInstructions.sh"%(os.environ["SCRATCH"],self.testName),"a+")
        self.CollectInstructions2File = open("%s/%s/plotInstructions.sh"%(os.environ["SCRATCH"],self.testName),"a+")


    def __WriteAlgorithmInfoForPlotting(self,TestID,AlgID)
        for param in self.AlgorithmList[TestID][0][AlgID]:
            ...write(..) echo "echo \"\${arg}\"" >> $SCRATCH/${testName}/plotInstructions.sh


    def __algorithmDispatch(self,TestID,AlgID,BinaryPath,scaleIndex,launchIndex,node,ppn,tpr):
        """
	"""
        # Set up the file string that will store the local benchmarking results
        fileString="DataFiles/%s_%dtest"%(self.AlgorithmList[TestID][0][AlgID].Tag,TestID)\
	          +"".join("_"+str(x) for x in self.AlgorithmList[TestID][0][AlgID].CurrentScaleParameters) + "%dlaunch_%dnodes_%dppn_%dtpr"%()
        # 'PreFile' requires NumNodes specification because in the 'Pre' stage, we want to keep the data for different node counts separate.
        PreFile="%s_%dtest"%(self.AlgorithmList[TestID][0][AlgID].Tag,TestID)\
               +"".join("_"+str(x) for x in self.AlgorithmList[TestID][0][AlgID].CurrentScaleParameters) + "%dlaunch_%dnodes_%dppn_%dtpr"%(....)
        PostFile="%s_%dtest"%(self.AlgorithmList[TestID][0][AlgID].Tag,TestID)\
               +"".join("_"+str(x) for x in self.AlgorithmList[TestID][0][AlgID].CurrentStartParameters) + "%dlaunch_%dppn_%dtpr"%(....)

        #UpdatePlotFile1="${tag1}_${scale}_${matrixDimMorig}_${matrixDimNorig}_${matrixDimKorig}_${cubeDimorig}"
        #UpdatePlotFile2="${tag1}_${scale}_${matrixDimMorig}__${matrixDimNorig}_${matrixDimKorig}${ppn}_${tpr}"

        # Plot instructions only need a single output per scaling study
        if (scaleIndex == 0):
            #WriteMethodDataForPlotting 0 ${UpdatePlotFile1} ${UpdatePlotFile2} ${tag1} ${PostFile} ${cubeDim} ${ppn} ${tpr}
            WriteAlgorithmInfoForPlotting(0,TestID,AlgID,launchIndex,ppn,tpr)	# Note that NumNodes is not included
            writePlotFileName ${PostFile} $SCRATCH/${testName}/plotInstructions.sh 1

        WriteAlgorithmInfoForCollectingStage1 ${tag1} ${PreFile} ${PreFile}_perf ${PreFile}_numerics $SCRATCH/${testName}/collectInstructionsStage1.sh
        WriteAlgorithmInfoForCollectingStage2 ${launchID} ${tag1} ${PreFile} ${PreFile}_perf ${PreFile}_numerics ${PostFile} ${PostFile}_perf ${PostFile}_numerics $SCRATCH/${testName}/collectInstructionsStage2.sh
        # Don't pass in 'cubeDim', because this is inferred based on the number of processes, as its just the cube root
        launchJobsPortal ${BinaryPath} ${tag1} ${fileString} ${curLaunchID} ${NumNodes} ${ppn} ${tpr} ${gemmORtrmm} ${algChoice} ${matrixDimM} ${matrixDimN} ${matrixDimK} ${numIterations} $SCRATCH/${testName}/${fileString}
        writePlotFileName ${fileString} $SCRATCH/${testName}/collectInstructionsStage1.sh 0


    def __portal(self,op,TestStartIndex,TestEndIndex,AlgIndex=0):
        """
        # Submit all scripts
        .. due to new format, need to watch out for repeats
        .. also note that different calls (3 so far) might want different things
        ..     writing scripts requires looking over every unique possible combination
        ..     launching scripts requires looking over every unique possible combination
        ..     calling the alg portal requires looking only over the test-specific node,ppn,tpr
        .. one solution here is to pass in a start,end testID parameters, and then use that as the outer-most loop
        .. and also to use a dictionary of tuples (node,ppn,tpr)
        """

        for LaunchIndex in range(1,NumLaunchesPerBinary+1):
            PortalDict = {}
            for TestIndex in range(TestStartIndex,TestEndIndex):
                curNumNodes=self.nodeMinList[TestIndex]
		scaleIndex=0
                while (curNumNodes <= self.nodeMaxList[TestIndex]):
                    curPPN=self.ppnMinList[TestIndex]
                    while (curPPN <= self.ppnMaxList[TestIndex]):
                        curTPR=self.tprMinList[TestIndex]
                        while (curTPR <= self.tprMaxList[TestIndex])
                            # Make sure we are in a suitable range
                            numPEsPerNode=curPPN*curTPR
                            if (minPEcountPerNode <= numPEsPerNode) and (maxPEcountPerNode >= numPEsPerNode):
                                add to PortalDict
                                if (op == 0):
                                    scriptName="%s/script_%s_round%s_launch%s_node%s_ppn%s_tpr%s.%s"%(self.testName,self.roundID,LaunchIndex,curNode,curPPN,curTPR,self.MachineType.BatchFileExtension)
                                    self.MachineType.queue(scriptName)
                                elif (op == 1):
                                    scriptName="%s/%s/script_%s_round%s_launch%s_node%s_ppn%s_tpr%s.%s"%(os.environ["SCRATCH"],self.fileID,self.roundID,LaunchIndex,curNode,curPPN,curTPR,self.MachineType.BatchFileExtension)
                                    scriptFile=open(scriptName,"a+")
                                    self.MachineType.script(scriptFile,self.testName,curNumNodes,curPPN,curTPR,numPEsPerNode,self.numHours,self.numMinutes,self.numSeconds)
                                elif (op == 2):
                                    algorithmDispatch(TestIndex,AlgIndex,BinaryPath,scaleIndex,LaunchIndex,curNumNodes,curPPN,curTPR):
                                .. what about checking in PortalDict??
                            curTPR=self.tprScaleOperatorList[TestIndex](curTPR,self.tprScaleFactorList[TestIndex])
                        curPPN=self.ppnScaleOperatorList[TestIndex](curPPN,self.ppnScaleFactorList[TestIndex])
                    scaleIndex=scaleIndex+1
                    curNumNodes=self.nodeScaleOperatorList[TestIndex](curNumNodes,self.nodeScaleFactorList[TestIndex])
		    if (op == 2):
		        self.AlgorithmList[TestIndex][0][AlgIndex].scale(TestIndex)

    def queue_submit(self):
        """
        """
	if (self.SubmitToQueue == 1):
          # Create directory to hold all binaries and then move them from ../Tests/testName/bin
          call("mkdir %s/%s/bin"%(os.environ["SCRATCH"],self.testName))
	  call("mv ../Tests/%s/bin/* %s/%s/bin"%(self.testName,os.environ["SCRATCH"],self.testName))
          portal(0,0,self.NumTests)

    def build(self):
        """
        """
        for lib in self.LibraryTypeList:
            os.environ["PROFTYPE"]="PERFORMANCE"
	    profType="P"
            # export SPECIAL_SCALA_ARG=REF
            lib.build()
            if (self.analyzeDecision1 == 1):
                profType="PC"
                os.environ["PROFTYPE"]="CRITTER"
                lib.build()

            if (self.analyzeDecision2 == 1):
                profType=profType+"T"
                export PROFTYPE=PROFILE
                os.environ["PROFTYPE"]="PROFILE"
                lib.build()

    def launch(self):
        """
        """
        os.environ["BINARYPATH"] = os.environ["SCRATCH"] + "/%s/bin/"%(self.testName)

        # collectData.sh will always be a single line, just a necessary intermediate step.
        ## echo "bash $SCRATCH/${testName}/collectInstructionsStage1.sh | bash PackageDataRemoteStage1.sh" > collectData.sh

        # Launch the generated script
        #bash $SCRATCH/${testName}.sh

        portal(1,0,self.NumTests)

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

        for TestIndex in range(0,numTests):
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

            for AlgIndex in range(len()):
                print("\nAlgorithm %s\n"%(self.AlgorithmList[TestIndex][0][AlgIndex].Tag))

                # Echo for SCAPLOT makefile generator
                binaryTag=self.AlgorithmList[...].Tag
                echo "echo \"\${binaryTag}\"" >> $SCRATCH/${testName}/collectInstructionsStage1.sh
                echo "echo \"\${binaryTag}\"" >> $SCRATCH/${testName}/collectInstructionsStage2.sh
                echo "echo \"\${binaryTag}\"" >> $SCRATCH/${testName}/plotInstructions.sh

                .. figure out what is going on with binaryPath and how it relates to experiments/machine/
                ..binaryPath=${BINARYPATH}\${binaryTag}_${machineName}
                if [ "${machineName}" == "PORTER" ];
                    binaryPath=\${binaryPath}_${mpiType}
                elif [ "${machineName}" == "BLUEWATERS" ];
                    # special case, only for CAMFS, not for bench_scalapack routines
                    if [ "\${binaryTag}" == "camfs_cacqr2" ] || [ "\${binaryTag}" == "camfs_cfr3d" ];
                        binaryPath=${BINARYPATH}\${binaryTag}_${machineName}_${GPU}

                portal(2,TestIndex,TestIndex,AlgIndex)
            echo "echo \"1\"" >> $SCRATCH/${testName}/collectInstructionsStage1.sh	# Signals end of the data files for this specific methodID
            echo "echo \"1\"" >> $SCRATCH/${testName}/collectInstructionsStage2.sh	# Signals end of the data files for this specific methodID
            echo "echo \"1\"" >> $SCRATCH/${testName}/plotInstructions.sh	# Signals end of the data files for this specific methodID

        queue_submit()

