import os
from subprocess import call
import datetime

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
        dateStr=datetime.datetime.now().strftime('%b-%d-%I%M%p-%G')
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
        os.environ["BINARYPATH"] = os.environ["SCRATCH"] + "/%s/bin"%(self.testName)

        call("mkdir %s/%s/"%(os.environ["SCRATCH"],self.testName),shell=True)
        call("mkdir %s/%s/DataFiles/"%(os.environ["SCRATCH"],self.testName),shell=True)

        self.PlotInstructionsFile = open("%s/%s/plotInstructions.sh"%(os.environ["SCRATCH"],self.testName),"a+")
        self.CollectInstructionsStage1File = open("%s/%s/plotInstructions.sh"%(os.environ["SCRATCH"],self.testName),"a+")
        self.CollectInstructionsStage2File = open("%s/%s/plotInstructions.sh"%(os.environ["SCRATCH"],self.testName),"a+")

    def __GetRangeCount(self,op,File,Start,End,Factor,Operator):
        """
	"""
	Curr = Start
	Count = 0
	while (Curr <= End):
	    if (op == 1):
	        File.write(str(Start)+"\n")
            Curr = Operator(Curr,Factor)
	    Count = Count + 1
	return Count

    def __WriteAlgorithmInfoForPlotting(self,TestID,AlgID)
        for param in self.AlgorithmList[TestID][0][AlgID]:
            self.PlotInstructionsFile.write(str(param)+"\n")

    def __WriteHeaderForCollection(self,File):
        """
        Writes the beginning of collectInstructionsStage1 and collectInstructionsStage2
	"""
        File.write("%s\n"%(self.testName))
        File.write("%s\n"%(self.testNameAllRounds))
        File.write("%s\n"%(self.MachineType.machineName))
        File.write("%d\n"%(self.numTests))


WriteMethodDataForCollectingStage1 () {
  local MethodTag=\${1}
  local FileNameBase=\${2}
  local FileName1=\${3}
  local FileName2=\${4}
  local WriteFile=\${5}

  echo "echo \"0\"" >> \${WriteFile}
  echo "echo \"\${MethodTag}\"" >> \${WriteFile}
  echo "echo \"\${FileName1}\"" >> \${WriteFile}
  if [ "\${MethodTag}" != "bscf" ];
  then
    echo "echo \"\${FileName2}\"" >> \${WriteFile}
  fi
  if [ "${profType}" == "PC" ] || [ "${profType}" == "PCT" ];
  then
    echo "echo \"\${FileNameBase}_critter\"" >> \${WriteFile}
  fi
  if [ "${profType}" == "PT" ] || [ "${profType}" == "PCT" ];
  then
    echo "echo \"\${FileNameBase}_timer\"" >> \${WriteFile}
  fi
}

WriteMethodDataForCollectingStage2 () {
  # Because 'Pre' (Stage1) collapses the NumLaunchesPerBinary, we do not want to overcount.
  local launchID=\${1}
  if [ \${launchID} -eq 1 ];
  then
    local MethodTag=\${2}
    local PreFileNameBase=\${3}
    local PreFileName1=\${4}
    local PreFileName2=\${5}
    local PostFileNameBase=\${6}
    local PostFileName1=\${7}
    local PostFileName2=\${8}
    local WriteFile=\${9}

    echo "echo \"0\"" >> \${WriteFile}
    echo "echo \"\${MethodTag}\"" >> \${WriteFile}
    echo "echo \"\${PostFileName1}\"" >> \${WriteFile}
    if [ "\${MethodTag}" != "bscf" ];
    then
      echo "echo \"\${PostFileName2}\"" >> \${WriteFile}
    fi
    if [ "${profType}" == "PC" ] || [ "${profType}" == "PCT" ];
    then
      echo "echo \"\${PostFileNameBase}_critter\"" >> \${WriteFile}
    fi
    if [ "${profType}" == "PT" ] || [ "${profType}" == "PCT" ];
    then
      echo "echo \"\${PostFileNameBase}_timer\"" >> \${WriteFile}
    fi
    echo "echo \"\${PreFileName1}\"" >> \${WriteFile}
    if [ "\${MethodTag}" != "bscf" ];
    then
      echo "echo \"\${PreFileName2}\"" >> \${WriteFile}
    fi
    if [ "${profType}" == "PC" ] || [ "${profType}" == "PCT" ];
    then
      echo "echo \"\${PreFileNameBase}_critter\"" >> \${WriteFile}
    fi
    if [ "${profType}" == "PT" ] || [ "${profType}" == "PCT" ];
    then
      echo "echo \"\${PreFileNameBase}_timer\"" >> \${WriteFile}
    fi
  fi
}


    def __writePlotFileName(self,..):
        # Performance runs will always run, so no reason for an if-statement here
        Prefix1=""
        Prefix2=""
        if [ "\${3}" == "1" ];
        then
            Prefix1="Raw/"
            Prefix2="Stats/"
        fi
        echo "echo \"\${Prefix1}\${1}_perf.txt\"" >> \${2}
        if [ "\${3}" == "1" ];
        then
            echo "echo \"\${Prefix2}\${1}_perf_stats.txt\"" >> \${2}
        fi
  
        echo "echo \"\${Prefix1}\${1}_numerics.txt\"" >> \${2}
        if [ "\${3}" == "1" ];
        then
            echo "echo \"\${Prefix2}\${1}_numerics_stats.txt\"" >> \${2}
        fi

        if [ "${profType}" == "PC" ] || [ "${profType}" == "PCT" ];
        then
            echo "echo \"\${1}_critter.txt\"" >> \${2}
            if [ "\${3}" == "1" ];
            then
                echo "echo \"\${1}_critter_breakdown.txt\"" >> \${2}
            fi
        fi
        if [ "${profType}" == "PT" ] || [ "${profType}" == "PCT" ];
        then
            echo "echo \"\${1}_timer.txt\"" >> \${2}
        fi


# Only for bsqr/rsqr -- only necessary for Performance now. Might want to use Critter later, but not Profiler
writePlotFileNameScalapackQR() {
  Prefix1=""
  Prefix2=""
  if [ "\${3}" == "1" ];
  then
    Prefix1="Raw/"
    Prefix2="Stats/"
  fi
  echo "echo \"\${Prefix1}\${1}_NoFormQ.txt\"" >> \${2}
  if [ "\${3}" == "1" ];
  then
    echo "echo \"\${Prefix2}\${1}_NoFormQ_stats.txt\"" >> \${2}
  fi
  echo "echo \"\${Prefix1}\${1}_FormQ.txt\"" >> \${2}
  if [ "\${3}" == "1" ];
  then
    echo "echo \"\${Prefix2}\${1}_FormQ_stats.txt\"" >> \${2}
  fi
}

# Only for bscf -- only necessary for Performance now. Might want to use Critter later, but not Profiler
writePlotFileNameScalapackCholesky() {
  Prefix1=""
  Prefix2=""
  if [ "\${3}" == "1" ];
  then
    Prefix1="Raw/"
    Prefix2="Stats/"
  fi
  echo "echo \"\${Prefix1}\${1}.txt\"" >> \${2}
  if [ "\${3}" == "1" ];
  then
    echo "echo \"\${Prefix2}\${1}_stats.txt\"" >> \${2}
  fi
}


# Functions that write the actual script, depending on machine
launchJobs () {
  local launchID=\${3}
  local numNodes=\${4}
  local ppn=\${5}
  local tpr=\${6}
  local numProcesses=\$((\${numNodes} * \${ppn}))
  local scriptName=$SCRATCH/${testName}/script_${fileID}id_${roundID}round_\${launchID}launchID_\${numNodes}nodes_\${ppn}ppn_\${tpr}tpr.${BatchFileExtension}
  writeTest \${numProcesses} \${ppn} \${tpr} \${scriptName} \${@:7:\$#}
}


launchJobsPortal () {
  # Launch performance job always.
  launchJobs \${@:2:6} \${1}_PERFORMANCE \${@:8:\${#}}

  # If analysis is turned on, launch Profiling job and Critter job.
  if [ "${profType}" == "PC" ] || [ "${profType}" == "PCT" ];
  then
    launchJobs \${@:2:6} \${1}_CRITTER \${@:8:\${#}}
  fi
  if [ "${profType}" == "PT" ] || [ "${profType}" == "PCT" ];
  then
    launchJobs \${@:2:6} \${1}_TIMER \${@:8:\${#}}
  fi
}




    def __algorithmDispatch(self,TestID,AlgID,BinaryPath,scaleIndex,launchIndex,node,ppn,tpr):
        """
	"""
        # Set up the file string that will store the local benchmarking results
        fileString="DataFiles/%s_%dtest"%(self.AlgorithmList[TestID][0][AlgID].Tag,TestID)\
	          +"".join("_"+str(x) for x in self.AlgorithmList[TestID][0][AlgID].CurrentScaleParameters) + "%dlaunch_%dnodes_%dppn_%dtpr"%(launchIndex,node,ppn,tpr)
        # 'PreFile' requires NumNodes specification because in the 'Pre' stage, we want to keep the data for different node counts separate.
        PreFile="%s_%dtest"%(self.AlgorithmList[TestID][0][AlgID].Tag,TestID)\
               +"".join("_"+str(x) for x in self.AlgorithmList[TestID][0][AlgID].CurrentScaleParameters) + "%dlaunch_%dnodes_%dppn_%dtpr"%(launchIndex,node,ppn,tpr)
        PostFile="%s_%dtest"%(self.AlgorithmList[TestID][0][AlgID].Tag,TestID)\
               +"".join("_"+str(x) for x in self.AlgorithmList[TestID][0][AlgID].CurrentStartParameters) + "%dlaunch_%dppn_%dtpr"%(launchIndex,node,ppn,tpr)

        #UpdatePlotFile1="${tag1}_${scale}_${matrixDimMorig}_${matrixDimNorig}_${matrixDimKorig}_${cubeDimorig}"
        #UpdatePlotFile2="${tag1}_${scale}_${matrixDimMorig}__${matrixDimNorig}_${matrixDimKorig}${ppn}_${tpr}"

	PrePath="%s/%s"%(os.environ["SCRATCH"],self.testName)
        # Plot instructions only need a single output per scaling study
        if (scaleIndex == 0):
            #WriteMethodDataForPlotting 0 ${UpdatePlotFile1} ${UpdatePlotFile2} ${tag1} ${PostFile} ${cubeDim} ${ppn} ${tpr}
            WriteAlgorithmInfoForPlotting(0,TestID,AlgID,launchIndex,ppn,tpr)	# Note that NumNodes is not included
            writePlotFileName ${PostFile} self.PlotInstructionsFile 1

        WriteAlgorithmInfoForCollectingStage1 ${tag1} ${PreFile} ${PreFile}_perf ${PreFile}_numerics self.CollectInstructionsStage1File
        WriteAlgorithmInfoForCollectingStage2 ${launchID} ${tag1} ${PreFile} ${PreFile}_perf ${PreFile}_numerics ${PostFile} ${PostFile}_perf ${PostFile}_numerics self.CollectInstructionsStage2File
        # Don't pass in 'cubeDim', because this is inferred based on the number of processes, as its just the cube root
        launchJobsPortal(BinaryPath,tag1,fileString,launchIndex,node,ppn,tpr,self.AlgorithmList[TestID][0][AlgID].CurrentScaleParameters,PrePath+"/%s"%(fileString))
        writePlotFileName fileString self.CollectInstructionsStage1File 0


    def __portal(self,op,TestStartIndex,TestEndIndex,AlgIndex=0):
        """
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
                                if (op == 0):
                                    TupleKey=(LaunchIndex,curNumNodes,curPPN,curTPR)
				    if not(TupleKey in PortalDict):
                                        scriptName="%s/script_%s_round%s_launch%s_node%s_ppn%s_tpr%s.%s"%(self.testName,self.roundID,LaunchIndex,curNode,curPPN,curTPR,self.MachineType.BatchFileExtension)
                                        self.MachineType.queue(scriptName)
				        PortalDict[TupleKey]=1
                                elif (op == 1):
                                    TupleKey=(LaunchIndex,curNumNodes,curPPN,curTPR)
				    if not(TupleKey in PortalDict):
                                        scriptName="%s/%s/script_%s_round%s_launch%s_node%s_ppn%s_tpr%s.%s"%(os.environ["SCRATCH"],self.fileID,self.roundID,LaunchIndex,curNode,curPPN,curTPR,self.MachineType.BatchFileExtension)
                                        scriptFile=open(scriptName,"a+")
                                        self.MachineType.script(scriptFile,self.testName,curNumNodes,curPPN,curTPR,numPEsPerNode,self.numHours,self.numMinutes,self.numSeconds)
				        PortalDict[TupleKey]=1
                                elif (op == 2):
                                    algorithmDispatch(TestIndex,AlgIndex,BinaryPath,scaleIndex,LaunchIndex,curNumNodes,curPPN,curTPR):
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
            # export SPECIAL_SCALA_ARG=REF
            lib.build()
            if (self.analyzeDecision1 == 1):
                os.environ["PROFTYPE"]="CRITTER"
                lib.build()

            if (self.analyzeDecision2 == 1):
                export PROFTYPE=PROFILE
                os.environ["PROFTYPE"]="PROFILE"
                lib.build()

    def launch(self):
        """
        """
        # collectData.sh will always be a single line, just a necessary intermediate step.
        ## echo "bash $SCRATCH/${testName}/collectInstructionsStage1.sh | bash PackageDataRemoteStage1.sh" > collectData.sh

        # Launch the generated script
        #bash $SCRATCH/${testName}.sh

        portal(1,0,self.NumTests)

        self.PlotInstructionsFile.write("1\n")
	self.PlotInstructionsFile.write("%s\n"%(self.testNameAllRounds))
	self.PlotInstructionsFile.write("%d\n"%(self.numTests))
	self.PlotInstructionsFile.write("%s\n"%(self.MachineType.machineName))
        WriteHeaderForCollection(self.CollectInstructionsStage1File)
        WriteHeaderForCollection(self.CollectInstructionsStage2File)

        for TestIndex in range(0,numTests):
            print("Test %d\n"%(i))

            self.PlotInstructionsFile.write("%s"%(self.AlgorithmList[TestIndex][1]))
            self.PlotInstructionsFile.write("%d"%(GetRangeCount(self.PlotInstructionsFile,0,self.nodeMinList[TestIndex],self.self.nodeMaxList[TestIndex],self.ppnScaleFactorList[TestIndex],self.ppnScaleOperatorList[TestIndex])))
            GetRangeCount(1,self.PlotInstructionsFile,self.nodeMinList[TestIndex],self.self.nodeMaxList[TestIndex],self.ppnScaleFactorList[TestIndex],self.ppnScaleOperatorList[TestIndex])

            # Figure out what PlotInstructions needs before fixing below
            #matrixDimM,matrixDimN,numIterations, comes from AlgList
            #echo "echo \"\${matrixDimM}\"" >> $SCRATCH/${testName}/plotInstructions.sh
            #echo "echo \"\${matrixDimN}\"" >> $SCRATCH/${testName}/plotInstructions.sh

            for AlgIndex in range(len()):
                print("\nAlgorithm %s\n"%(self.AlgorithmList[TestIndex][0][AlgIndex].Tag))

                # Echo for SCAPLOT makefile generator
                binaryTag=self.AlgorithmList[...].Tag
                self.CollectInstructionsStage1File.write("%s\n"%(self.AlgorithmList[TestIndex][0][AlgIndex].Tag))
                self.CollectInstructionsStage2File.write("%s\n"%(self.AlgorithmList[TestIndex][0][AlgIndex].Tag))
                self.PlotInstructionsFile.write("%s\n"%(self.AlgorithmList[TestIndex][0][AlgIndex].Tag))

                binaryPath="%s/%s"%(os.environ["BINARYPATH"],self.AlgorithmList[TestIndex][0][AlgIndex].Tag)
                # Below: special case that will hopefully be replaced soon
                if (self.MachineType.IsAccelerated())
                    binaryPath=binaryPath + "_GPU"

                portal(2,TestIndex,TestIndex,AlgIndex)
            self.CollectInstructionsStage1File.write("1\n")
            self.CollectInstructionsStage2File.write("1\n")
            self.PlotInstructionsFile.write("1\n")
        queue_submit()

