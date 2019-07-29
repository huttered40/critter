import os

class bench(object):
    """
    """
    def __init__(self,\
                 CritterPath,\
		 MachineType,\
		 LibraryTypeList,\
		 fileID,\
		 roundID,\
		 minNumNodes,\
		 maxNumNodes,\
		 nodeScaleFactor,\
		 ppnScaleFactor,\
		 tprScaleFactor,\
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
		 ppnMinList,\
		 ppnMaxList,\
		 tprMinList,\
		 tprMaxList,\
		 SubmitToQueue,\
		 AlgorithmList):
       """
       CritterPath - specify full path to Critter repository
                   - do not include a '/' after 'critter'

       MachineType - machine type specified in ../experiments/machines/

       LibraryTypeList - list of types specified in ../experiments/libraries/

       fileID - base name of the directory inside which all data/scripts will be stored
              - will appear inside the SCRATCH directory

       roundID - set to '1' unless performing piecewise testing (launching same job separately) to enhance performance reproducibility

       minNumNodes - minimum number of nodes needed for any one test

       maxNumNodes - maximum number of nodes needed for any one test

       nodeScaleFactor - scaling factor to apply to the number of nodes

       ppnScaleFactor - scaling factor to apply to the number of MPI processes per node (ppn)

       tprScaleFactor - scaling factor to apply to the number of threads per MPI rank (tpr)

       NumLaunchesPerBinary - set to '1' unless performing testing to enhance performance reproducibility (by launching same jobs multiple times)
                            - different from 'roundID' because the former did not launch at the same time, but waited via a separate launch of bench.sh

       numTests - number of scaling studies
                - for example, weak scaling and strong scaling, even across the same variants, constitute separate tests

       numHours,numMinutes,numSeconds - ....

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

       ppnMinList,ppnMaxList,tprMinList,tprMaxList - Lists with an element for each node count
                                                   - specify the min and max ppn,tpr at each node count
                                                   - 'ppn' stands for 'MPI processes per node'
                                                   - 'tpr' stands for 'threads-per-MPI-rank'

       SubmitToQueue - '1' to submit jobs to queue, '0' to not submit to queue

       AlgorithmList - list of lists of lists of inputs, must be of length 'numTests'
                          - each inner list holds the inputs for the corresponding algorithm in AlgorithmList
       """
       self.CritterPath = CritterPath
       self.MachineType = MachineType
       self.LibraryTypeList = LibraryTypeList
       self.fileID = fileID
       self.roundID = roundID
       self.minNumNodes = minNumNodes
       self.maxNumNodes = maxNumNodes
       self.nodeScaleFactor = nodeScaleFactor
       self.ppnScaleFactor = ppnScaleFactor
       self.tprScaleFactor = tprScaleFactor
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
       self.ppnMinList = ppnMinList
       self.ppnMaxList = ppnMaxList
       self.tprMinList = tprMinList
       self.tprMaxList = tprMaxList
       self.SubmitToQueue = SubmitToQueue
       self.AlgorithmList = AlgorithmList

    def launch(self):
        """
        """
        self.MachineType.set()
