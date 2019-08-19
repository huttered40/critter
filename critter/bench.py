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
        call("mkdir %s/%s/bin"%(os.environ["SCRATCH"],self.testName),shell=True)

        self.PlotInstructionsFile = open("%s/%s/plotInstructions.txt"%(os.environ["SCRATCH"],self.testName),"a+")
        self.CollectInstructionsFile = open("%s/%s/collectInstructions.txt"%(os.environ["SCRATCH"],self.testName),"a+")

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

    def WriteAlgorithmInfoForPlotting(self,AlgParameters,launchID,ppn,tpr):
        for param in AlgParameters:
            self.PlotInstructionsFile.write(str(param)+"\n")
        self.PlotInstructionsFile.write(str(launchID)+"\n")
	self.PlotInstructionsFile.write(str(ppn)+"\n")
	self.PlotInstructionsFile.write(str(tpr)+"\n")

    def WriteHeaderForCollection(self,File):
        """
        Writes the beginning of collectInstructionsStage1 and collectInstructionsStage2
	"""
        File.write("%s\n"%(self.testName))
        File.write("%s\n"%(self.testNameAllRounds))
        File.write("%s\n"%(self.MachineType.MachineName))
        File.write("%d\n"%(self.numTests))

    def WriteAlgInfoForCollecting(self,launchID,TestID,AlgTag,PreFile,PostFile):
        if (launchID == 1):
            self.CollectInstructionsFile.write("0\n")
            self.CollectInstructionsFile.write("%s\n"%(AlgTag))
            FileExtensions=self.TestList[TestID][2]
            self.CollectInstructionsFile.write("%d\n"%(len(FileExtensions)))

            # Allow for any number of user-defined tests
	    for i in range(len(FileExtensions)):
                self.CollectInstructionsFile.write("%s_%s\n"%(PreFile,FileExtensions[i][0]))
                self.CollectInstructionsFile.write("%s_%s\n"%(PostFile,FileExtensions[i][0]))

    def writePlotFileName(self,TestID,DataFile,WriteFile,Tag):
        # Performance runs will always run, so no reason for an if-statement here
        Prefix1=""
        Prefix2=""
        if (Tag == 1):
            Prefix1="Raw"
            Prefix2="Stats"

        FileExtensions=self.TestList[TestID][2]
        # Allow for any number of user-defined tests
	for i in range(len(FileExtensions)):
            WriteFile.write("%s/%s_%s.txt\n"%(Prefix1,DataFile,FileExtensions[i][0]))
            if (Tag==1):
                WriteFile.write("%s/%s_%s_stats.txt\n"%(Prefix2,DataFile,FileExtensions[i][0]))

    # Functions that write the actual script, depending on machine
    def launchJobs(self,BinaryPath,launchIndex,TestID,AlgID,node,ppn,tpr,AlgParameters,fileString):
        """
	"""
	FileExtensions=self.TestList[TestID][2]
        numProcesses=node*ppn
        scriptName="%s/%s/script_%s_round%s_launch%s_node%s_ppn%s_tpr%s.%s"%(os.environ["SCRATCH"],self.testName,self.fileID,self.roundID,launchIndex,node,ppn,tpr,self.MachineType.BatchFileExtension)

        # Allow for any number of user-defined tests
        MethodString = BinaryPath+"".join(" "+str(x) for x in AlgParameters);
	for i in range(len(FileExtensions)):
            MethodString = MethodString + " %s_%s"%(fileString,FileExtensions[i][0])
        scriptFile=open(scriptName,"a+")
	self.MachineType.write_test(scriptFile,numProcesses,ppn,tpr,MethodString)
        scriptFile.close()

    def algorithmDispatch(self,TestID,AlgParameters,AlgID,BinaryPath,scaleIndex,launchID,node,ppn,tpr):
        """
	"""
        # Set up the file string that will store the local benchmarking results
        BaseString="%s_%dtest"%(self.TestList[TestID][0][AlgID].Tag,TestID)\
                  +"".join("_"+str(x) for x in AlgParameters) + "_%dlaunch_%dppn_%dtpr"%(launchID,ppn,tpr)
        PostFile=BaseString
        # 'PreFile' requires NumNodes specification because in the 'Pre' stage, we want to keep the data for different node counts separate.
        PreFile=BaseString+"_%dnodes"%(node)
        fileString="DataFiles/"+PreFile

        #UpdatePlotFile1="${tag1}_${scale}_${matrixDimMorig}_${matrixDimNorig}_${matrixDimKorig}_${cubeDimorig}"
        #UpdatePlotFile2="${tag1}_${scale}_${matrixDimMorig}__${matrixDimNorig}_${matrixDimKorig}${ppn}_${tpr}"

	PrePath="%s/%s"%(os.environ["SCRATCH"],self.testName)
        # Plot instructions only need a single output per scaling study
        if (scaleIndex == 0):
            # look at position of the UpdatePlotFile* files WriteMethodDataForPlotting 0 ${UpdatePlotFile1} ${UpdatePlotFile2} ${tag1} ${PostFile} ${cubeDim} ${ppn} ${tpr}
            self.WriteAlgorithmInfoForPlotting(AlgParameters,launchID,ppn,tpr)	# Note that NumNodes is not included
            self.writePlotFileName(TestID,PostFile,self.PlotInstructionsFile,1)
            pass

        self.WriteAlgInfoForCollecting(launchID,TestID,self.TestList[TestID][0][AlgID].Tag,PreFile,PostFile)
        self.launchJobs(BinaryPath,launchID,TestID,AlgID,node,ppn,tpr,AlgParameters,PrePath+"/%s"%(fileString))
        #self.writePlotFileName(TestID,PostFile,self.CollectInstructionsFile,0)


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
                                        scriptName="%s/script_%s_round%s_launch%s_node%s_ppn%s_tpr%s.%s"%(self.testName,self.fileID,self.roundID,LaunchIndex,curNumNodes,curPPN,curTPR,self.MachineType.BatchFileExtension)
                                        self.MachineType.queue(scriptName)
				        PortalDict[TupleKey]=1
                                elif (op == 1):
                                    TupleKey=(LaunchIndex,curNumNodes,curPPN,curTPR)
				    if not(TupleKey in PortalDict):
                                        scriptName="%s/%s/script_%s_round%s_launch%s_node%s_ppn%s_tpr%s.%s"%(os.environ["SCRATCH"],self.testName,self.fileID,self.roundID,LaunchIndex,curNumNodes,curPPN,curTPR,self.MachineType.BatchFileExtension)
                                        scriptFile=open(scriptName,"a+")
                                        self.MachineType.script(scriptFile,self.testName,curNumNodes,curPPN,curTPR,numPEsPerNode,self.numHours,self.numMinutes,self.numSeconds)
                                        scriptFile.close()
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
	# Create directory to hold all binaries and then move them from ../Tests/testName/bin
	call("mv %s/Tests/%s/bin/* %s/%s/bin"%(self.CritterPath,self.testName,os.environ["SCRATCH"],self.testName),shell=True)
        self.portal(0,0,self.numTests)

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
            self.CollectInstructionsFile.write("%s\n"%(BinaryTag))
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
        self.portal(1,0,self.numTests)

        self.PlotInstructionsFile.write("1\n")
	self.PlotInstructionsFile.write("%s\n"%(self.testNameAllRounds))
	self.PlotInstructionsFile.write("%d\n"%(self.numTests))
	self.PlotInstructionsFile.write("%s\n"%(self.MachineType.MachineName))
        self.WriteHeaderForCollection(self.CollectInstructionsFile)

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
            self.CollectInstructionsFile.write("1\n")
            self.PlotInstructionsFile.write("1\n")

    def launch(self):
        self.queue_submit()

