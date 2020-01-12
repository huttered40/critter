import os,sys
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
		 CritterBreakdownInfo,\
		 CritterVizInfo,\
		 fileID,\
		 roundID,\
		 NumLaunchesPerBinary,\
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

        CritterBreakdownInfo - string that matches the critical_path_breakdown bitset global variable in src/critter.h

        CritterVizInfo - list that gives info as how to how use scaplot visualization tool to plot and predict data
                       - [0] - Use CritterViz (generate micro benchmarks)(2), Output to file(1), or output to standard output(0)
                       - [1] - Min message size for benchmarking the MPI collectives
                       - [2] - Max message size for benchmarking the MPI collectives
                       - [3] - Jump factor for message size for benchmarking the MPI collectives
                       - [4] - Min subcommunicator size for benchmarking the MPI collectives
                       - [5] - Jump factor for subcommunicator size for benchmarking the MPI collectives

        fileID - base name of the directory inside which all data/scripts will be stored
               - will appear inside the SCRATCH directory
               - specified as a string

        roundID - set to '1' unless performing piecewise testing (launching same job separately) to enhance performance reproducibility

        NumLaunchesPerBinary - set to '1' unless performing testing to enhance performance reproducibility (by launching same jobs multiple times)
                             - different from 'roundID' because the former did not launch at the same time, but waited via a separate launch of bench.sh

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
        assert(len(CritterBreakdownInfo)==6)
        self.CritterPath = CritterPath
        self.MachineType = MachineType
        self.CritterBreakdownInfo = CritterBreakdownInfo
        self.UseCritterViz = CritterVizInfo[0]
        if (self.UseCritterViz>0):
            self.MinMessageSize = CritterVizInfo[1]
            self.MaxMessageSize = CritterVizInfo[2]
            self.MessageJumpSize = CritterVizInfo[3]
            self.MinProcessSize = CritterVizInfo[4]
            # Note: no MaxProcessSize because this is dependent on the (node,ppn,tpr)
            self.ProcessJumpSize = CritterVizInfo[5]
            self.num_iter = CritterVizInfo[6]
        self.fileID = fileID
        self.roundID = roundID
        self.NumLaunchesPerBinary = NumLaunchesPerBinary
        self.numTests = len(TestList)
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

        if (self.UseCritterViz>0):
	    # First check if the user correctly set the corresponding environment variable. (It cannot be set inside this script!)
            if ("CRITTER_VIZ" not in os.environ.keys()):#os.environ["CRITTER_VIZ"]==""):
                print("User must set the environment variable `CRITTER_VIZ`=`ON`")
                sys.exit(0)
            if (self.UseCritterViz>1):
                # Build the micro benchmarks and move all binaries from ../tests/testName/bin
                call("cd %s; make test"%(self.CritterPath),shell=True)
        # I think these directories serve mainly as a intermediate place to put the binaries
	#   before being moved to SCRATCH
        call("mkdir %s/tests/%s"%(self.CritterPath,self.testName),shell=True)
        call("mkdir %s/tests/%s/bin"%(self.CritterPath,self.testName),shell=True)
        # Copy all binaries in bin/ into test folder's bin
        call("cp %s/bin/* %s/tests/%s/bin/"%(self.CritterPath,self.CritterPath,self.testName),shell=True)

        self.MachineType.set()
        os.environ["BINARYPATH"] = os.environ["SCRATCH"] + "/%s/bin"%(self.testName)

        call("mkdir %s/%s/"%(os.environ["SCRATCH"],self.testName),shell=True)
        call("mkdir %s/%s/data/"%(os.environ["SCRATCH"],self.testName),shell=True)
        call("mkdir %s/%s/bin"%(os.environ["SCRATCH"],self.testName),shell=True)

        self.PlotInstructionsFile = open("%s/%s/plotInstructions.txt"%(os.environ["SCRATCH"],self.testName),"a+")
        self.CollectInstructionsFile = open("%s/%s/collectInstructions.txt"%(os.environ["SCRATCH"],self.testName),"a+")
        self.BenchInstructionsFile = open("%s/%s/benchInstructions.txt"%(os.environ["SCRATCH"],self.testName),"a+")

    def stringify(self,x):
        """
        fixes a file-error when negative integers are casted into string characters that cannot be specified when opening files
        expects an integer (for now)
        """
        if (x<0):
            return "-"+str((-1)*x)
        else:
            return str(x)

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
                if (self.TestList[TestIndex][AlgIndex].SpecialFunc(AlgParameterList,[curNumNodes,curPPN,curTPR])):
		    totalScaleCount+=1
            curNumNodes=self.nodeScaleOperatorList[TestIndex](curNumNodes,self.nodeScaleFactorList[TestIndex])
            self.TestList[TestIndex][AlgIndex].scale(AlgParameterList,scaleIndex)
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
        File.write("%s\n"%(self.CritterBreakdownInfo))

    def WriteAlgInfoForPlotting(self,AlgParameters,launchID,ppn,tpr):
        if (launchID==1):
            self.PlotInstructionsFile.write(str(len(AlgParameters))+"\n")
            for param in AlgParameters:
                self.PlotInstructionsFile.write(self.stringify(param)+"\n")
	    self.PlotInstructionsFile.write(str(ppn)+"\n")
	    self.PlotInstructionsFile.write(str(tpr)+"\n")

    def WriteAlgInfoForCollecting(self,launchID,File,TestID,AlgID,AlgTag,PreFile,PostFile):
        if (launchID == 1):
            File.write("0\n")
            File.write("%s\n"%(AlgTag))
            File.write("%d\n"%(4+len(self.TestList[TestID][AlgID].InputParameterStartRange)+2))	# '+5' numPEs,testID,launchID,ppn,tpr,(node?)
            File.write("%s+critter.txt\n"%(PreFile))
            File.write("%s+critter.txt\n"%(PostFile))

    def WriteAlgInfoForBenchmarking(self,launchID,File,AlgTag,PreFile,PostFile):
        if (launchID == 1):
            File.write("%s\n"%(AlgTag))
            File.write("%s.txt\n"%(PreFile))
            File.write("%s.txt\n"%(PostFile))

    # Functions that write the actual script, depending on machine
    def launchJobs(self,BinaryPath,launchIndex,TestID,AlgID,node,ppn,tpr,AlgParameters,fileString):
        """
	"""
        numProcesses=node*ppn
        scriptName="%s/%s/script_%s_round%s_launch%s_node%s_ppn%s_tpr%s.%s"%(os.environ["SCRATCH"],self.testName,self.fileID,self.roundID,launchIndex,node,ppn,tpr,self.MachineType.BatchFileExtension)

        # Allow for any number of user-defined tests
        MethodString = BinaryPath+"".join(" "+str(x) for x in AlgParameters)#+" %d %d"%(ppn,tpr);
        scriptFile=open(scriptName,"a+")
        scriptFile.write("export CRITTER_VIZ_FILE=%s+critter\n"%(fileString))
	self.MachineType.write_test(scriptFile,numProcesses,ppn,tpr,MethodString)
        scriptFile.close()

    def write_benchmark(self,scriptFile,BinaryPath,launchID,nodes,ppn,tpr):
        """
	"""
	# For now, we want to test 5 MPI routines (see below)
	for tag in [("test_send_recv",0),("test_sendrecv",0),("test_bcast",0),("test_reduce",0),("test_allreduce",0),("test_allgather",1),("test_sendrecv_replace",0)]:
            BinaryPath="%s/%s"%(os.environ["BINARYPATH"],tag[0])
            if (tag[1] == 0):	# routines whose message size does not depend on process count
                msg_size=self.MinMessageSize
                while (msg_size <= self.MaxMessageSize):
                    pcount=self.MinProcessSize;
                    while (pcount<=(nodes*ppn)):
                        AlgParameters=[msg_size,pcount,self.num_iter]
                        # Set up the file string that will store the local benchmarking results
                        BaseString1="%s"%(tag[0])\
                            +"".join("+"+str(x) for x in AlgParameters) + "+%d+%d+%d"%(launchID,ppn,tpr)
                        BaseString2="%s"%(tag[0])\
                            + "+%d+%d+%d"%(launchID,ppn,tpr)
                        PostFile=BaseString2
                        # 'PreFile' requires NumNodes specification because in the 'Pre' stage, we want to keep the data for different node counts separate.
                        PreFile=BaseString1+"+%d"%(nodes)
	                PrePath="%s/%s"%(os.environ["SCRATCH"],self.testName)
                        fileString=PrePath+"/data/"+PreFile
                        MethodString = BinaryPath+"".join(" "+str(x) for x in AlgParameters)#+" %d %d"%(ppn,tpr);
                        self.WriteAlgInfoForBenchmarking(launchID,self.BenchInstructionsFile,tag[0],PreFile,PostFile)
                        scriptFile.write("export CRITTER_VIZ_FILE=%s\n"%(fileString))
	                self.MachineType.write_test(scriptFile,nodes*ppn,ppn,tpr,MethodString)
                        pcount*=self.ProcessJumpSize;
                    msg_size*=self.MessageJumpSize
            else:
                pcount=self.MinProcessSize;
                while (pcount<=(nodes*ppn)):
                    msg_size=self.MinMessageSize*pcount
                    while (msg_size <= self.MaxMessageSize):
                        AlgParameters=[msg_size,pcount,self.num_iter]
                        # Set up the file string that will store the local benchmarking results
                        BaseString1="%s"%(tag[0])\
                            +"".join("+"+str(x) for x in AlgParameters) + "+%d+%d+%d"%(launchID,ppn,tpr)
                        BaseString2="%s"%(tag[0])\
                            +"".join("+"+str(x) for x in self.SaveAlgParameters) + "+%d+%d+%d"%(launchID,ppn,tpr)
                        PostFile=BaseString2
                        # 'PreFile' requires NumNodes specification because in the 'Pre' stage, we want to keep the data for different node counts separate.
                        PreFile=BaseString1+"+%d"%(nodes)
	                PrePath="%s/%s"%(os.environ["SCRATCH"],self.testName)
                        fileString=PrePath+"/data/"+PreFile
                        MethodString = BinaryPath+"".join(" "+str(x) for x in AlgParameters)#+" %d %d"%(ppn,tpr);
                        self.WriteAlgInfoForBenchmarking(launchID,self.BenchInstructionsFile,tag[0],PreFile,PostFile)
                        scriptFile.write("export CRITTER_VIZ_FILE=%s\n"%(fileString))
	                self.MachineType.write_test(scriptFile,nodes*ppn,ppn,tpr,MethodString)
                        msg_size*=self.MessageJumpSize
                    pcount*=self.ProcessJumpSize;

    def algorithmDispatch(self,TestID,AlgParameters,AlgID,BinaryPath,IsFirstNode,scaleIndex,launchID,node,ppn,tpr):
        """
	"""
        # Set up the file string that will store the local benchmarking results
        BaseString1="%s+%d"%(self.TestList[TestID][AlgID].Tag,TestID)\
                  +"".join("+"+self.stringify(x) for x in AlgParameters) + "+%d+%d+%d"%(launchID,ppn,tpr)
        BaseString2="%s+%d"%(self.TestList[TestID][AlgID].Tag,TestID)\
                  +"".join("+"+self.stringify(x) for x in self.SaveAlgParameters) + "+%d+%d+%d"%(launchID,ppn,tpr)
        PostFile=BaseString2
        # 'PreFile' requires NumNodes specification because in the 'Pre' stage, we want to keep the data for different node counts separate.
        PreFile=BaseString1+"+%d"%(node)
        fileString="data/"+PreFile

	PrePath="%s/%s"%(os.environ["SCRATCH"],self.testName)
        # Plot instructions only need a single output per scaling study
        if (IsFirstNode):
            self.SavePostFile = PostFile
            # look at position of the UpdatePlotFile* files WriteMethodDataForPlotting 0 ${UpdatePlotFile1} ${UpdatePlotFile2} ${tag1} ${PostFile} ${cubeDim} ${ppn} ${tpr}
            self.WriteAlgInfoForCollecting(launchID,self.PlotInstructionsFile,TestID,AlgID,self.TestList[TestID][AlgID].Tag,PreFile,PostFile)
            self.WriteAlgInfoForPlotting(self.SaveAlgParameters,launchID,ppn,tpr)	# Note that NumNodes is not included
        self.WriteAlgInfoForCollecting(launchID,self.CollectInstructionsFile,TestID,AlgID,self.TestList[TestID][AlgID].Tag,PreFile,PostFile)
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
                        if (op == 2):
                            ValidNodeList.append([])
                            ValidProcessList.append([])
                        curNumNodes=self.GetNodeListOffset(TestIndex,0)
                        PrevAlgParameterList = list(AlgParameterList)
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
                                    if (self.TestList[TestIndex][AlgIndex].SpecialFunc(AlgParameterList,[curNumNodes,curPPN,curTPR])):
                                        TupleKey=(LaunchIndex,curNumNodes,curPPN,curTPR)
				        if not(TupleKey in self.PortalDict):
                                            scriptName="%s/%s/script_%s_round%s_launch%s_node%s_ppn%s_tpr%s.%s"%(os.environ["SCRATCH"],self.testName,self.fileID,self.roundID,LaunchIndex,curNumNodes,curPPN,curTPR,self.MachineType.BatchFileExtension)
                                            scriptFile=open(scriptName,"a+")
                                            self.MachineType.script(scriptFile,self.testName,curNumNodes,curPPN,curTPR,numPEsPerNode,self.numHours,self.numMinutes,self.numSeconds)
                                            if (self.UseCritterViz==2):
                                                self.write_benchmark(scriptFile,BinaryPath,LaunchIndex,curNumNodes,curPPN,curTPR)
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
                                        PrevAlgParameterList = list(AlgParameterList)
                                    else:
                                        print("Variant with wrong params - ", AlgParameterList,[curNumNodes,curPPN,curTPR])
                                        AlgParameterList = list(PrevAlgParameterList)
		            if (op == 2):
                                self.TestList[TestIndex][AlgIndex].scale(AlgParameterList,scaleIndex)
                            scaleIndex=scaleIndex+1
                            curNumNodes=self.nodeScaleOperatorList[TestIndex](curNumNodes,self.nodeScaleFactorList[TestIndex])
		        if (op == 2):
                            if (ValidNodeList[-1] == []):
                                del ValidNodeList[-1]
                                del ValidProcessList[-1]
                        curTPR=self.tprScaleOperatorList[TestIndex](curTPR,self.tprScaleFactorList[TestIndex])
                    curPPN=self.ppnScaleOperatorList[TestIndex](curPPN,self.ppnScaleFactorList[TestIndex])

    def queue_submit(self):
        """
        """
        call("mv %s/tests/%s/bin/* %s/%s/bin"%(self.CritterPath,self.testName,os.environ["SCRATCH"],self.testName),shell=True)
        self.portal(0,0,1)

    def cycle(self,TestIndex,AlgIndex,VariantIndex,ParameterIndex,AlgParameterList,ValidNodeList,ValidProcessList):
        """
	"""
        # base case at the last level -- this is when we are sure a valid parameter combination exists
	if (ParameterIndex == self.TestList[TestIndex][AlgIndex].NumParameters):
            # Echo for SCAPLOT makefile generator
            BinaryTag=self.TestList[TestIndex][AlgIndex].Tag

            BinaryPath="%s/%s"%(os.environ["BINARYPATH"],BinaryTag)
            print("\n    Variant %d"%(VariantIndex))
            self.portal(2,TestIndex,TestIndex+1,list(AlgParameterList),AlgIndex,BinaryPath,ValidNodeList,ValidProcessList)
	    return VariantIndex+1

        IsValid=1
	while (IsValid):
	    if (AlgParameterList[ParameterIndex] == self.TestList[TestIndex][AlgIndex].InputParameterEndRange[ParameterIndex]):
	        IsValid=0
            VariantIndex = self.cycle(TestIndex,AlgIndex,VariantIndex,ParameterIndex+1,list(AlgParameterList),ValidNodeList,ValidProcessList)
            if (IsValid):
	        AlgParameterList[ParameterIndex] = self.TestList[TestIndex][AlgIndex].InputParameterScaleOperator[ParameterIndex](\
                    AlgParameterList[ParameterIndex],self.TestList[TestIndex][AlgIndex].InputParameterScaleFactor[ParameterIndex])
	    else:
	        return VariantIndex

    def generate(self):
        """
        """

        self.WriteHeader(self.PlotInstructionsFile)
        self.WriteHeader(self.CollectInstructionsFile)
        self.PlotInstructionsFile.write(str(self.MachineType.PeakNodePerformance)+"\n")
        self.PlotInstructionsFile.write(str(self.MachineType.PeakNetworkInjectionRate)+"\n")

        self.PortalDict = {}
        for TestIndex in range(0,self.numTests):
            print("\nTest %d"%(TestIndex))

            # These two lists below will have length equal to the number of valid algorithm variants
            ValidNodeList=[]
            ValidProcessList=[]
            #self.PortalDict = {}
            for AlgIndex in range(len(self.TestList[TestIndex])):
                print("\n  Algorithm %s"%(self.TestList[TestIndex][AlgIndex].Tag))
                VariantIndex=0
		AlgParameterList=list(self.TestList[TestIndex][AlgIndex].InputParameterStartRange)
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

