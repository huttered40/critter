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

        ppnMinList,ppnMaxList - min/max number of processes-per-node for each test

        tprMinList,tprMaxList - min/max number of threads-per-process for each test

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

        # These list is necessary for special tracking
        self.SaveAlgParameters = []
        self.SaveJobStrDict = {}
        self.NodeCountDict = {}
        self.ProcessCountDict = {}
        self.PortalDict = {}

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
	        File.write(str(Curr)+"\n")
            Curr = Operator(Curr,Factor)
	    Count = Count + 1
	return Count

    def GetNodeListOffset(self,TestIndex,StartNodeIndex):
        cur=self.nodeMinList[TestIndex]
        counter=0
        while (counter < StartNodeIndex):
            cur = self.nodeScaleOperatorList[TestIndex](cur,self.nodeScaleFactorList[TestIndex])
            counter+=1
        return cur

    def GetTotalValidNodes(self,TestIndex,AlgIndex,AlgParameterList,LaunchIndex,curPPN,curTPR,StartNodeIndex):
        """
        """
        totalScaleCount=0
        scaleIndex=0
        curNumNodes=self.GetNodeListOffset(TestIndex,StartNodeIndex)
        while (curNumNodes <= self.nodeMaxList[TestIndex]):
            # Make sure we are in a suitable range
	    numPEsPerNode=curPPN*curTPR
            if (self.minPEcountPerNode <= numPEsPerNode) and (self.maxPEcountPerNode >= numPEsPerNode):
                if (self.TestList[TestIndex][0][AlgIndex].SpecialFunc(AlgParameterList,[curNumNodes,curPPN,curTPR])):
		    totalScaleCount+=1
            curNumNodes=self.nodeScaleOperatorList[TestIndex](curNumNodes,self.nodeScaleFactorList[TestIndex])
            self.TestList[TestIndex][0][AlgIndex].scale(AlgParameterList,scaleIndex)
            scaleIndex+=1
        return totalScaleCount

    def WriteHeader(self,File):
        """
        Writes the beginning of collectInstructionsStage1 and collectInstructionsStage2
	"""
        File.write("%s\n"%(self.testName))
        File.write("%s\n"%(self.testNameAllRounds))
        File.write("%s\n"%(self.MachineType.MachineName))
        File.write("%d\n"%(self.numTests))

    def WriteAlgorithmInfoForPlotting(self,AlgParameters,launchID,ppn,tpr):
        if (launchID==1):
            self.PlotInstructionsFile.write(str(len(AlgParameters))+"\n")
            for param in AlgParameters:
                self.PlotInstructionsFile.write(str(param)+"\n")
	    self.PlotInstructionsFile.write(str(ppn)+"\n")
	    self.PlotInstructionsFile.write(str(tpr)+"\n")

    def WriteAlgInfoForCollecting(self,launchID,File,TestID,AlgID,AlgTag,PreFile,PostFile):
        if (launchID == 1):
            File.write("0\n")
            File.write("%s\n"%(AlgTag))
            FileExtensions=self.TestList[TestID][2]
            File.write("%d\n"%(len(FileExtensions)))
            File.write("%d\n"%(1+2*(len(self.TestList[TestID][0][AlgID].InputParameterStartRange)+2)))	# '+2' from ppn,tpr

            # Allow for any number of user-defined tests
	    for i in range(len(FileExtensions)):
                File.write("%s_%s.txt\n"%(PreFile,FileExtensions[i][0]))
                File.write("%s_%s.txt\n"%(PostFile,FileExtensions[i][0]))

    # Functions that write the actual script, depending on machine
    def launchJobs(self,BinaryPath,launchIndex,TestID,AlgID,node,ppn,tpr,AlgParameters,fileString):
        """
	"""
	FileExtensions=self.TestList[TestID][2]
        numProcesses=node*ppn
        scriptName="%s/%s/script_%s_round%s_launch%s_node%s_ppn%s_tpr%s.%s"%(os.environ["SCRATCH"],self.testName,self.fileID,self.roundID,launchIndex,node,ppn,tpr,self.MachineType.BatchFileExtension)

        # Allow for any number of user-defined tests
        MethodString = BinaryPath+"".join(" "+str(x) for x in AlgParameters)+" %d %d"%(ppn,tpr);
	for i in range(len(FileExtensions)):
            MethodString = MethodString + " %s_%s"%(fileString,FileExtensions[i][0])
        scriptFile=open(scriptName,"a+")
	self.MachineType.write_test(scriptFile,numProcesses,ppn,tpr,MethodString)
        scriptFile.close()

    def algorithmDispatch(self,TestID,AlgParameters,AlgID,BinaryPath,IsFirstNode,scaleIndex,launchID,node,ppn,tpr):
        """
	"""
        # Set up the file string that will store the local benchmarking results
        BaseString1="%s_%dtest"%(self.TestList[TestID][0][AlgID].Tag,TestID)\
                  +"".join("_"+str(x) for x in AlgParameters) + "_%dlaunch_%dppn_%dtpr"%(launchID,ppn,tpr)
        BaseString2="%s_%dtest"%(self.TestList[TestID][0][AlgID].Tag,TestID)\
                  +"".join("_"+str(x) for x in self.SaveAlgParameters) + "_%dlaunch_%dppn_%dtpr"%(launchID,ppn,tpr)
        PostFile=BaseString2
        # 'PreFile' requires NumNodes specification because in the 'Pre' stage, we want to keep the data for different node counts separate.
        PreFile=BaseString1+"_%dnodes"%(node)
        fileString="DataFiles/"+PreFile

	PrePath="%s/%s"%(os.environ["SCRATCH"],self.testName)
        # Plot instructions only need a single output per scaling study
        if (IsFirstNode):
            self.SavePostFile = PostFile
            # look at position of the UpdatePlotFile* files WriteMethodDataForPlotting 0 ${UpdatePlotFile1} ${UpdatePlotFile2} ${tag1} ${PostFile} ${cubeDim} ${ppn} ${tpr}
            self.WriteAlgInfoForCollecting(launchID,self.PlotInstructionsFile,TestID,AlgID,self.TestList[TestID][0][AlgID].Tag,PreFile,PostFile)
            self.WriteAlgorithmInfoForPlotting(self.SaveAlgParameters,launchID,ppn,tpr)	# Note that NumNodes is not included
        self.WriteAlgInfoForCollecting(launchID,self.CollectInstructionsFile,TestID,AlgID,self.TestList[TestID][0][AlgID].Tag,PreFile,PostFile)
        self.launchJobs(BinaryPath,launchID,TestID,AlgID,node,ppn,tpr,AlgParameters,PrePath+"/%s"%(fileString))


    def portal(self,op,TestStartIndex,TestEndIndex,AlgParameterList=[],AlgIndex=0,BinaryPath=0,ValidNodeList=[],ValidProcessList=[]):
        """
        Note that 'portal' exploits the fact that as we scale, the range of PPN/TPR counts will not change.
        However, this does not mean that the PPN/TPR counts are the same for each node.
        The same algorithm variant can start at different node counts and thus reach an entirely different collection of node counts with the same NumNodes multiplier
        """
	SaveAlgParameterList = list(AlgParameterList)
        for LaunchIndex in range(1,self.NumLaunchesPerBinary+1):
            for TestIndex in range(TestStartIndex,TestEndIndex):
                curPPN = self.ppnMinList[TestIndex]
                while (curPPN <= self.ppnMaxList[TestIndex]):
                    curTPR = self.tprMinList[TestIndex]
                    while (curTPR <= self.tprMaxList[TestIndex]):
		        scaleIndex=0
                        IsFirstNode=True
			# Must reset 'AlgParameterList' each time to avoid corrupting its elements as they are modified across nodes
	                AlgParameterList = list(SaveAlgParameterList)
                        # Two lines below assume that this algorithm variant will have at least one valid node
                        ValidNodeList.append([])
                        ValidProcessList.append([])
                        curNumNodes=self.GetNodeListOffset(TestIndex,0)
                        while (curNumNodes <= self.nodeMaxList[TestIndex]):
                            # Make sure we are in a suitable range
	                    numPEsPerNode=curPPN*curTPR
                            if (self.minPEcountPerNode <= numPEsPerNode) and (self.maxPEcountPerNode >= numPEsPerNode):
                                if (op == 0):
                                    TupleKey=(LaunchIndex,curNumNodes,curPPN,curTPR)
	            	            if (TupleKey in self.SaveJobStrDict):
                                        scriptName="%s/script_%s_round%s_launch%s_node%s_ppn%s_tpr%s.%s"%(self.testName,self.fileID,self.roundID,LaunchIndex,curNumNodes,curPPN,curTPR,self.MachineType.BatchFileExtension)
                                        self.MachineType.queue(scriptName)
				        self.PortalDict[TupleKey]=1
                                elif (op == 2):
                                    # Save special variables if at 1st node count
                                    if (scaleIndex == 0):
                                        self.SaveAlgParameters = list(AlgParameterList)
                                    if (self.TestList[TestIndex][0][AlgIndex].SpecialFunc(AlgParameterList,[curNumNodes,curPPN,curTPR])):
                                        TupleKey=(LaunchIndex,curNumNodes,curPPN,curTPR)
				        if not(TupleKey in self.PortalDict):
                                            scriptName="%s/%s/script_%s_round%s_launch%s_node%s_ppn%s_tpr%s.%s"%(os.environ["SCRATCH"],self.testName,self.fileID,self.roundID,LaunchIndex,curNumNodes,curPPN,curTPR,self.MachineType.BatchFileExtension)
                                            scriptFile=open(scriptName,"a+")
                                            self.MachineType.script(scriptFile,self.testName,curNumNodes,curPPN,curTPR,numPEsPerNode,self.numHours,self.numMinutes,self.numSeconds)
                                            scriptFile.close()
				            self.PortalDict[TupleKey]=1
                                            self.SaveJobStrDict[TupleKey]=1
			                self.algorithmDispatch(TestIndex,AlgParameterList,AlgIndex,BinaryPath,IsFirstNode,scaleIndex,LaunchIndex,curNumNodes,curPPN,curTPR)
                                        IsFirstNode=False
                                        # Save the node/process counts supporting valid variants
                                        self.NodeCountDict[curNumNodes]=1
                                        self.ProcessCountDict[curNumNodes*curPPN]=1
                                        ValidNodeList[-1].append(curNumNodes)
                                        ValidProcessList[-1].append(curNumNodes*curPPN)
                                    else:
                                        print("Variant with wrong params - ", AlgParameterList,[curNumNodes,curPPN,curTPR])
		            if (op == 2):
                                self.TestList[TestIndex][0][AlgIndex].scale(AlgParameterList,scaleIndex)
                            scaleIndex=scaleIndex+1
                            curNumNodes=self.nodeScaleOperatorList[TestIndex](curNumNodes,self.nodeScaleFactorList[TestIndex])
                        curTPR=self.tprScaleOperatorList[TestIndex](curTPR,self.tprScaleFactorList[TestIndex])
                    curPPN=self.ppnScaleOperatorList[TestIndex](curPPN,self.ppnScaleFactorList[TestIndex])

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
            lib.build(self.CritterPath,self.testName)

    def cycle(self,TestIndex,AlgIndex,VariantIndex,ParameterIndex,AlgParameterList,ValidNodeList,ValidProcessList):
        """
	"""
        # base case at the last level -- this is when we are sure a valid parameter combination exists
	if (ParameterIndex == self.TestList[TestIndex][0][AlgIndex].NumParameters):
            # Echo for SCAPLOT makefile generator
            BinaryTag=self.TestList[TestIndex][0][AlgIndex].Tag

            BinaryPath="%s/%s"%(os.environ["BINARYPATH"],BinaryTag)
            print("\n    Variant %d"%(VariantIndex))
            self.portal(2,TestIndex,TestIndex+1,list(AlgParameterList),AlgIndex,BinaryPath,ValidNodeList,ValidProcessList)
	    return VariantIndex+1

        IsValid=1
	while (IsValid):
	    if (AlgParameterList[ParameterIndex] == self.TestList[TestIndex][0][AlgIndex].InputParameterEndRange[ParameterIndex]):
	        IsValid=0
            VariantIndex = self.cycle(TestIndex,AlgIndex,VariantIndex,ParameterIndex+1,list(AlgParameterList),ValidNodeList,ValidProcessList)
            if (IsValid):
	        AlgParameterList[ParameterIndex] = self.TestList[TestIndex][0][AlgIndex].InputParameterScaleOperator[ParameterIndex](\
                    AlgParameterList[ParameterIndex],self.TestList[TestIndex][0][AlgIndex].InputParameterScaleFactor[ParameterIndex])
	    else:
	        return VariantIndex

    def generate(self):
        """
        """
        self.WriteHeader(self.PlotInstructionsFile)
        self.WriteHeader(self.CollectInstructionsFile)
        self.PlotInstructionsFile.write(str(self.MachineType.PeakNodePerformance)+"\n")
        self.PlotInstructionsFile.write(str(self.MachineType.PeakNetworkInjectionRate)+"\n")

        for TestIndex in range(0,self.numTests):
            print("\nTest %d"%(TestIndex))

            self.PlotInstructionsFile.write("%s\n"%(self.TestList[TestIndex][1]))
            self.PlotInstructionsFile.write("%d\n"%(len(self.TestList[TestIndex][2])))
            if (self.TestList[TestIndex][2][0][0] == "critter"):
                self.PlotInstructionsFile.write("1\n")
                self.CollectInstructionsFile.write("1\n")
                NonCritterIndex=1
            else:
                self.PlotInstructionsFile.write("0\n")
                self.CollectInstructionsFile.write("0\n")
                NonCritterIndex=0
            self.PlotInstructionsFile.write("%d\n"%(len(self.TestList[TestIndex][2][NonCritterIndex][1])))
            for ColumnHeader in self.TestList[TestIndex][2][NonCritterIndex][1]:
                self.PlotInstructionsFile.write("%s\n"%(ColumnHeader))

            # These two lists below will have length equal to the number of valid algorithm variants
            ValidNodeList=[]
            ValidProcessList=[]
            self.PortalDict = {}
            for AlgIndex in range(len(self.TestList[TestIndex][0])):
                print("\n  Algorithm %s"%(self.TestList[TestIndex][0][AlgIndex].Tag))
                VariantIndex=0
		AlgParameterList=list(self.TestList[TestIndex][0][AlgIndex].InputParameterStartRange)
		self.cycle(TestIndex,AlgIndex,VariantIndex,0,AlgParameterList,ValidNodeList,ValidProcessList)

            # Signify end of test info
            self.CollectInstructionsFile.write("1\n")
            self.PlotInstructionsFile.write("1\n")
            # Detail number and names of the valid node counts and process counts
            count=0
	    self.PlotInstructionsFile.write("%d\n"%(len(self.NodeCountDict.keys())))
            for key in sorted(self.NodeCountDict.keys()):
                self.NodeCountDict[key] = count
                self.PlotInstructionsFile.write("%d\n"%(key))
                count+=1
	    self.PlotInstructionsFile.write("%d\n"%(len(self.ProcessCountDict.keys())))
            count=0
            for key in sorted(self.ProcessCountDict.keys()):
                self.ProcessCountDict[key] = count
                self.PlotInstructionsFile.write("%d\n"%(key))
                count+=1
            # Detail the nodeScaleIndices and processScaleIndices for each algorithm variant
            self.PlotInstructionsFile.write("%d\n"%(len(ValidNodeList)))
            for i in range(len(ValidNodeList)):
                self.PlotInstructionsFile.write("%d\n"%(len(ValidNodeList[i])))
                for j in range(len(ValidNodeList[i])):
                    self.PlotInstructionsFile.write("%d\n"%(self.NodeCountDict[ValidNodeList[i][j]]))
            self.PlotInstructionsFile.write("%d\n"%(len(ValidProcessList)))
            for i in range(len(ValidProcessList)):
                self.PlotInstructionsFile.write("%d\n"%(len(ValidProcessList[i])))
                for j in range(len(ValidProcessList[i])):
                    self.PlotInstructionsFile.write("%d\n"%(self.ProcessCountDict[ValidProcessList[i][j]]))
            self.NodeCountDict.clear()
            self.ProcessCountDict.clear()

    def launch(self):
        self.queue_submit()
        self.CollectInstructionsFile.close()
        self.PlotInstructionsFile.close()

