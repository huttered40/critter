import os
from subprocess import call
import datetime

# See if this works with the lambda specified in instructions2.py
# I don't really want to include all possible libraries that the user might use in the lambda
import math

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
		 TestList):
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

        analyzeDecision1 - profile using Critter

        analyzeDecision2 - profile using TAU
                         - not currently supported

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

        TestList - list of lists
	              - each inner list holds:
		          - a list of 'algorithm' instances to be tested against
			  - a string detailing the kind of scaling study
			  - a list holding a list of strings detailing the specific launch, where each string is accompanied by a tuple of column names used for plotting,
			     each listed in the same data file
			        - strings must be in order corresponding to how they are called in source code.
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
        self.analyzeDecision1 = analyzeDecision1
        self.analyzeDecision2 = analyzeDecision2
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
        self.TestList = TestList
        dateStr=datetime.datetime.now().strftime('%b-%d-%I%M%p-%G')
        self.testName="%s_%s_%s_round%d"%(fileID,dateStr,self.MachineType.MachineName,roundID)
        self.testNameAllRounds="%s_%s"%(fileID,self.MachineType.MachineName)

        # TODO: create conditional check to see if Tests/ exists. Its not tracked by git.
	call("mkdir %s/Tests"%(self.CritterPath),shell=True)

        # I think these directories serve mainly as a intermediate place to put the binaries
	#   before being moved to SCRATCH
        call("mkdir %s/Tests/%s"%(self.CritterPath,self.testName),shell=True)
        call("mkdir %s/Tests/%s/bin"%(self.CritterPath,self.testName),shell=True)

        self.MachineType.set()
        os.environ["BINARYPATH"] = os.environ["SCRATCH"] + "/%s/bin"%(self.testName)

        call("mkdir %s/%s/"%(os.environ["SCRATCH"],self.testName),shell=True)
        call("mkdir %s/%s/DataFiles/"%(os.environ["SCRATCH"],self.testName),shell=True)

        self.PlotInstructionsFile = open("%s/%s/plotInstructions.sh"%(os.environ["SCRATCH"],self.testName),"a+")
        self.CollectInstructionsStage1File = open("%s/%s/collectInstructionsStage1.sh"%(os.environ["SCRATCH"],self.testName),"a+")
        self.CollectInstructionsStage2File = open("%s/%s/collectInstructionsStage2.sh"%(os.environ["SCRATCH"],self.testName),"a+")

    def GetRangeCount(self,op,File,Start,End,Factor,Operator):
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

    def WriteAlgorithmInfoForPlotting(self,TestID,AlgID):
        for param in self.TestList[TestID][0][AlgID]:
            self.PlotInstructionsFile.write(str(param)+"\n")

    def WriteHeaderForCollection(self,File):
        """
        Writes the beginning of collectInstructionsStage1 and collectInstructionsStage2
	"""
        File.write("%s\n"%(self.testName))
        File.write("%s\n"%(self.testNameAllRounds))
        File.write("%s\n"%(self.MachineType.MachineName))
        File.write("%d\n"%(self.numTests))

    """
    def __WriteMethodDataForCollectingStage1(self):
        MethodTag=\${1}
        FileNameBase=\${2}
        FileName1=\${3}
        FileName2=\${4}
        WriteFile=\${5}

        echo "echo \"0\"" >> \${WriteFile}
        echo "echo \"\${MethodTag}\"" >> \${WriteFile}
        echo "echo \"\${FileName1}\"" >> \${WriteFile}
        if [ "\${MethodTag}" != "bscf" ];
            echo "echo \"\${FileName2}\"" >> \${WriteFile}
        if [ "${profType}" == "PC" ] || [ "${profType}" == "PCT" ];
            echo "echo \"\${FileNameBase}_critter\"" >> \${WriteFile}
        if [ "${profType}" == "PT" ] || [ "${profType}" == "PCT" ];
            echo "echo \"\${FileNameBase}_timer\"" >> \${WriteFile}

    def __WriteMethodDataForCollectingStage2(self):
        # Because 'Pre' (Stage1) collapses the NumLaunchesPerBinary, we do not want to overcount.
        launchID=\${1}
        if [ \${launchID} -eq 1 ];
            MethodTag=\${2}
            PreFileNameBase=\${3}
            PreFileName1=\${4}
            PreFileName2=\${5}
            PostFileNameBase=\${6}
            PostFileName1=\${7}
            PostFileName2=\${8}
            WriteFile=\${9}

            echo "echo \"0\"" >> \${WriteFile}
            echo "echo \"\${MethodTag}\"" >> \${WriteFile}
            echo "echo \"\${PostFileName1}\"" >> \${WriteFile}
            if [ "\${MethodTag}" != "bscf" ];
                echo "echo \"\${PostFileName2}\"" >> \${WriteFile}
            if [ "${profType}" == "PC" ] || [ "${profType}" == "PCT" ];
                echo "echo \"\${PostFileNameBase}_critter\"" >> \${WriteFile}
            if [ "${profType}" == "PT" ] || [ "${profType}" == "PCT" ];
                echo "echo \"\${PostFileNameBase}_timer\"" >> \${WriteFile}
            echo "echo \"\${PreFileName1}\"" >> \${WriteFile}
            if [ "\${MethodTag}" != "bscf" ];
                echo "echo \"\${PreFileName2}\"" >> \${WriteFile}
            if [ "${profType}" == "PC" ] || [ "${profType}" == "PCT" ];
                echo "echo \"\${PreFileNameBase}_critter\"" >> \${WriteFile}
            if [ "${profType}" == "PT" ] || [ "${profType}" == "PCT" ];
                echo "echo \"\${PreFileNameBase}_timer\"" >> \${WriteFile}


    def __writePlotFileName(self,..):
        # Performance runs will always run, so no reason for an if-statement here
        Prefix1=""
        Prefix2=""
        if [ "\${3}" == "1" ];
            Prefix1="Raw/"
            Prefix2="Stats/"
        echo "echo \"\${Prefix1}\${1}_perf.txt\"" >> \${2}
        if [ "\${3}" == "1" ];
            echo "echo \"\${Prefix2}\${1}_perf_stats.txt\"" >> \${2}
  
        echo "echo \"\${Prefix1}\${1}_numerics.txt\"" >> \${2}
        if [ "\${3}" == "1" ];
            echo "echo \"\${Prefix2}\${1}_numerics_stats.txt\"" >> \${2}

        if [ "${profType}" == "PC" ] || [ "${profType}" == "PCT" ];
            echo "echo \"\${1}_critter.txt\"" >> \${2}
            if [ "\${3}" == "1" ];
                echo "echo \"\${1}_critter_breakdown.txt\"" >> \${2}
        if [ "${profType}" == "PT" ] || [ "${profType}" == "PCT" ];
            echo "echo \"\${1}_timer.txt\"" >> \${2}


    # Only for bsqr/rsqr -- only necessary for Performance now. Might want to use Critter later, but not Profiler
    def __writePlotFileNameScalapackQR(self):
	Prefix1=""
        Prefix2=""
        if [ "\${3}" == "1" ];
            Prefix1="Raw/"
            Prefix2="Stats/"
        echo "echo \"\${Prefix1}\${1}_NoFormQ.txt\"" >> \${2}
        if [ "\${3}" == "1" ];
            echo "echo \"\${Prefix2}\${1}_NoFormQ_stats.txt\"" >> \${2}
        echo "echo \"\${Prefix1}\${1}_FormQ.txt\"" >> \${2}
        if [ "\${3}" == "1" ];
            echo "echo \"\${Prefix2}\${1}_FormQ_stats.txt\"" >> \${2}


    # Only for bscf -- only necessary for Performance now. Might want to use Critter later, but not Profiler
    def __writePlotFileNameScalapackCholesky(self):
	Prefix1=""
        Prefix2=""
        if [ "\${3}" == "1" ];
            Prefix1="Raw/"
            Prefix2="Stats/"
        echo "echo \"\${Prefix1}\${1}.txt\"" >> \${2}
        if [ "\${3}" == "1" ];
            echo "echo \"\${Prefix2}\${1}_stats.txt\"" >> \${2}
    """

    # Functions that write the actual script, depending on machine
    def launchJobs(self,BinaryPath,launchIndex,node,ppn,tpr,AlgorithmInfo,AlgParameters,fileString):
        """
	"""
        numProcesses=node*ppn
        scriptName="%s/%s/script_%s_round%s_launch%s_node%s_ppn%s_tpr%s.%s"%(os.environ["SCRATCH"],self.testName,self.fileID,self.roundID,launchIndex,node,ppn,tpr,self.MachineType.BatchFileExtension)

        MethodString = BinaryPath+"".join(" "+str(x) for x in AlgParameters)+" %s"%(fileString)
	self.MachineType.writeTest(numProcesses,ppn,tpr,MethodString)

    def algorithmDispatch(self,TestID,AlgParameters,AlgID,BinaryPath,scaleIndex,launchIndex,node,ppn,tpr):
        """
	"""
        # Set up the file string that will store the local benchmarking results
        fileString="DataFiles/%s_%dtest"%(self.TestList[TestID][0][AlgID].Tag,TestID)\
	          +"".join("_"+str(x) for x in AlgParameters) + "_%dlaunch_%dnodes_%dppn_%dtpr"%(launchIndex,node,ppn,tpr)
        # 'PreFile' requires NumNodes specification because in the 'Pre' stage, we want to keep the data for different node counts separate.
        #PreFile="%s_%dtest"%(self.TestList[TestID][0][AlgID].Tag,TestID)\
        #       +"".join("_"+str(x) for x in self.TestList[TestID][0][AlgID].CurrentScaleParameters) + "%dlaunch_%dnodes_%dppn_%dtpr"%(launchIndex,node,ppn,tpr)
        #PostFile="%s_%dtest"%(self.TestList[TestID][0][AlgID].Tag,TestID)\
        #       +"".join("_"+str(x) for x in self.TestList[TestID][0][AlgID].CurrentStartParameters) + "%dlaunch_%dppn_%dtpr"%(launchIndex,node,ppn,tpr)

        #UpdatePlotFile1="${tag1}_${scale}_${matrixDimMorig}_${matrixDimNorig}_${matrixDimKorig}_${cubeDimorig}"
        #UpdatePlotFile2="${tag1}_${scale}_${matrixDimMorig}__${matrixDimNorig}_${matrixDimKorig}${ppn}_${tpr}"

	PrePath="%s/%s"%(os.environ["SCRATCH"],self.testName)
        # Plot instructions only need a single output per scaling study
        if (scaleIndex == 0):
            # look at position of the UpdatePlotFile* files WriteMethodDataForPlotting 0 ${UpdatePlotFile1} ${UpdatePlotFile2} ${tag1} ${PostFile} ${cubeDim} ${ppn} ${tpr}
            #WriteAlgorithmInfoForPlotting(0,TestID,AlgID,launchIndex,ppn,tpr)	# Note that NumNodes is not included
            #writePlotFileName ${PostFile} self.PlotInstructionsFile 1
            pass

        #WriteAlgorithmInfoForCollectingStage1 ${tag1} ${PreFile} ${PreFile}_perf ${PreFile}_numerics self.CollectInstructionsStage1File
        #WriteAlgorithmInfoForCollectingStage2 ${launchID} ${tag1} ${PreFile} ${PreFile}_perf ${PreFile}_numerics ${PostFile} ${PostFile}_perf ${PostFile}_numerics self.CollectInstructionsStage2File
        # Don't pass in 'cubeDim', because this is inferred based on the number of processes, as its just the cube root
        self.launchJobs(BinaryPath,launchIndex,node,ppn,tpr,self.TestList[TestID][0][AlgID],AlgParameters,PrePath+"/%s"%(fileString))
        #writePlotFileName fileString self.CollectInstructionsStage1File 0


    def portal(self,op,TestStartIndex,TestEndIndex,AlgParameterList=[],AlgIndex=0,BinaryPath=0):
        """
        """
        for LaunchIndex in range(1,self.NumLaunchesPerBinary+1):
            PortalDict = {}
            for TestIndex in range(TestStartIndex,TestEndIndex):
                curNumNodes=self.nodeMinList[TestIndex]
		scaleIndex=0
                while (curNumNodes <= self.nodeMaxList[TestIndex]):
                    curPPN=self.ppnMinList[TestIndex][scaleIndex]
                    while (curPPN <= self.ppnMaxList[TestIndex][scaleIndex]):
                        curTPR=self.tprMinList[TestIndex][scaleIndex]
                        while (curTPR <= self.tprMaxList[TestIndex][scaleIndex]):
                            # Make sure we are in a suitable range
			    numPEsPerNode=curPPN*curTPR
                            if (self.minPEcountPerNode <= numPEsPerNode) and (self.maxPEcountPerNode >= numPEsPerNode):
                                if (op == 0):
                                    TupleKey=(LaunchIndex,curNumNodes,curPPN,curTPR)
				    if not(TupleKey in PortalDict):
                                        scriptName="%s/script_%s_round%s_launch%s_node%s_ppn%s_tpr%s.%s"%(self.testName,self.fileID,self.roundID,LaunchIndex,curNode,curPPN,curTPR,self.MachineType.BatchFileExtension)
                                        self.MachineType.queue(scriptName)
				        PortalDict[TupleKey]=1
                                elif (op == 1):
                                    TupleKey=(LaunchIndex,curNumNodes,curPPN,curTPR)
				    if not(TupleKey in PortalDict):
                                        scriptName="%s/%s/script_%s_round%s_launch%s_node%s_ppn%s_tpr%s.%s"%(os.environ["SCRATCH"],self.testName,self.fileID,self.roundID,LaunchIndex,curNumNodes,curPPN,curTPR,self.MachineType.BatchFileExtension)
                                        scriptFile=open(scriptName,"a+")
                                        self.MachineType.script(scriptFile,self.testName,curNumNodes,curPPN,curTPR,numPEsPerNode,self.numHours,self.numMinutes,self.numSeconds)
				        PortalDict[TupleKey]=1
                                elif (op == 2):
                                    if (self.TestList[TestIndex][0][AlgIndex].SpecialFunc(AlgParameterList,[LaunchIndex,curNumNodes,curPPN,curTPR])):
				        self.algorithmDispatch(TestIndex,AlgParameterList,AlgIndex,BinaryPath,scaleIndex,LaunchIndex,curNumNodes,curPPN,curTPR)
                            curTPR=self.tprScaleOperatorList[TestIndex](curTPR,self.tprScaleFactorList[TestIndex])
                        curPPN=self.ppnScaleOperatorList[TestIndex](curPPN,self.ppnScaleFactorList[TestIndex])
                    curNumNodes=self.nodeScaleOperatorList[TestIndex](curNumNodes,self.nodeScaleFactorList[TestIndex])
		    if (op == 2):
		        self.TestList[TestIndex][0][AlgIndex].scale(AlgParameterList,scaleIndex)
                    scaleIndex=scaleIndex+1

    def queue_submit(self):
        """
        """
	if (self.SubmitToQueue == 1):
	    # Create directory to hold all binaries and then move them from ../Tests/testName/bin
            call("mkdir %s/%s/bin"%(os.environ["SCRATCH"],self.testName),shell=True)
	    call("mv ../Tests/%s/bin/* %s/%s/bin"%(self.testName,os.environ["SCRATCH"],self.testName),shell=True)
            self.portal(0,0,self.NumTests)

    def build(self):
        """
        """
        for lib in self.LibraryTypeList:
            # export SPECIAL_SCALA_ARG=REF
            lib.build(self.CritterPath,self.testName)

    def cycle(self,TestIndex,AlgIndex,VariantIndex,ParameterIndex,AlgParameterList):
        """
	"""
        # base case at the last level -- this is when we are sure a valid parameter combination exists
	if (ParameterIndex == self.TestList[TestIndex][0][AlgIndex].NumParameters):
            # Echo for SCAPLOT makefile generator
            BinaryTag=self.TestList[TestIndex][0][AlgIndex].Tag
            self.CollectInstructionsStage1File.write("%s\n"%(BinaryTag))
            self.CollectInstructionsStage2File.write("%s\n"%(BinaryTag))
            self.PlotInstructionsFile.write("%s\n"%(BinaryTag))

            BinaryPath="%s/%s"%(os.environ["BINARYPATH"],BinaryTag)
            # Below: special case that will hopefully be replaced soon
            if (self.MachineType.IsAccelerated()):
                BinaryPath=BinaryPath + "_GPU"

            print("\n    Variant %d"%(VariantIndex))
            self.portal(2,TestIndex,TestIndex+1,list(AlgParameterList),AlgIndex,BinaryPath)
	    return VariantIndex+1

        IsValid=1
	while (IsValid):
	    if (AlgParameterList[ParameterIndex] == self.TestList[TestIndex][0][AlgIndex].InputParameterEndRange[ParameterIndex]):
	        IsValid=0
            VariantIndex = self.cycle(TestIndex,AlgIndex,VariantIndex,ParameterIndex+1,list(AlgParameterList))
            if (IsValid):
	        AlgParameterList[ParameterIndex] = self.TestList[TestIndex][0][AlgIndex].InputParameterScaleOperator[ParameterIndex](\
                    AlgParameterList[ParameterIndex],self.TestList[TestIndex][0][AlgIndex].InputParameterScaleFactor[ParameterIndex])
	    else:
	        return VariantIndex

    def generate(self):
        """
        """
        # collectData.sh will always be a single line, just a necessary intermediate step.
        ## echo "bash $SCRATCH/${testName}/collectInstructionsStage1.sh | bash PackageDataRemoteStage1.sh" > collectData.sh

        # Launch the generated script
        #bash $SCRATCH/${testName}.sh

        self.portal(1,0,self.numTests)

        self.PlotInstructionsFile.write("1\n")
	self.PlotInstructionsFile.write("%s\n"%(self.testNameAllRounds))
	self.PlotInstructionsFile.write("%d\n"%(self.numTests))
	self.PlotInstructionsFile.write("%s\n"%(self.MachineType.MachineName))
        self.WriteHeaderForCollection(self.CollectInstructionsStage1File)
        self.WriteHeaderForCollection(self.CollectInstructionsStage2File)

        for TestIndex in range(0,self.numTests):
            print("\nTest %d"%(TestIndex))

            self.PlotInstructionsFile.write("%s"%(self.TestList[TestIndex][1]))
            NodeCount = self.GetRangeCount(0,self.PlotInstructionsFile,self.nodeMinList[TestIndex],self.nodeMaxList[TestIndex],self.ppnScaleFactorList[TestIndex],self.ppnScaleOperatorList[TestIndex])
	    self.PlotInstructionsFile.write("%d"%(NodeCount))
            self.GetRangeCount(1,self.PlotInstructionsFile,self.nodeMinList[TestIndex],self.nodeMaxList[TestIndex],self.ppnScaleFactorList[TestIndex],self.ppnScaleOperatorList[TestIndex])

            # Figure out what PlotInstructions needs before fixing below
            #matrixDimM,matrixDimN,numIterations, comes from AlgList
            #echo "echo \"\${matrixDimM}\"" >> $SCRATCH/${testName}/plotInstructions.sh
            #echo "echo \"\${matrixDimN}\"" >> $SCRATCH/${testName}/plotInstructions.sh

            for AlgIndex in range(len(self.TestList[TestIndex][0])):
                print("\n  Algorithm %s"%(self.TestList[TestIndex][0][AlgIndex].Tag))
                VariantIndex=0
		AlgParameterList=list(self.TestList[TestIndex][0][AlgIndex].InputParameterStartRange)
		self.cycle(TestIndex,AlgIndex,VariantIndex,0,AlgParameterList)
            self.CollectInstructionsStage1File.write("1\n")
            self.CollectInstructionsStage2File.write("1\n")
            self.PlotInstructionsFile.write("1\n")

    def launch(self):
        self.queue_submit()

