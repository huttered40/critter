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
                   - specified as a string

       MachineType - machine type specified in ../experiments/machines/

       LibraryTypeList - list of types specified in ../experiments/libraries/

       fileID - base name of the directory inside which all data/scripts will be stored
              - will appear inside the SCRATCH directory
              - specified as a string

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

machineName=""
accelType=""
testAccel_NoAccel=""
if [ "$(hostname |grep "porter")" != "" ];
  machineName=PORTER
elif [ "$(hostname |grep "stampede2")" != "" ];
  machineName=STAMPEDE2
elif [ "$(hostname |grep "h2o")" != "" ];
  accelType="n"
  testAccel_NoAccel="n"
  machineName=BLUEWATERS
  export GPU=NOGPU

if [ "${mpiType}" == "mpi" ];
  export MPITYPE=MPI_TYPE
elif [ "${mpiType}" == "ampi" ];
  export MPITYPE=AMPI_TYPE

dateStr=$(date +%Y-%m-%d-%H_%M_%S)
testName=${fileID}_${dateStr}_${machineName}_round${roundID}
testNameAllRounds=${fileID}_${machineName}	# Name of the corresponding directory in CAMFS_data. Allows for appending multiple runs
mkdir ../Tests/${testName}
mkdir ../Tests/${testName}/bin

if [ ${dataType} == 0 ];
  export DATATYPE=FLOAT_TYPE
elif [ ${dataType} == 1 ];
  export DATATYPE=DOUBLE_TYPE
elif [ ${dataType} == 2 ];
  export DATATYPE=COMPLEX_FLOAT_TYPE
elif [ ${dataType} == 3 ];
  export DATATYPE=COMPLEX_DOUBLE_TYPE
if [ ${intType} == 0 ];
  export INTTYPE=INT_TYPE
elif [ ${intType} == 1 ];
  export INTTYPE=INT64_T_TYPE


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
bash $SCRATCH/${testName}.sh

# Submit job scripts to queue - note that Porter doesn't need to do this
if [ "${SubmitToQueue}" == "1" ];
then
  mkdir $SCRATCH/${testName}/bin
  mv ../Tests/${testName}/* $SCRATCH/${testName}/bin
  cd $SCRATCH

  # Submit all scripts
  curLaunchID=1
  while [ ${curLaunchID} -le ${NumLaunchesPerBinary} ];
  do
    curNumNodes=${minNumNodes}
    listIndex=0
    while [ ${curNumNodes} -le ${maxNumNodes} ];
    do
      curPPN=${ppnMinList[${listIndex}]}
      ppnMax=${ppnMaxList[${listIndex}]}
      while [ ${curPPN} -le ${ppnMax} ];
      do
        curTPR=${tprMinList[${listIndex}]}
        tprMax=${tprMaxList[${listIndex}]}
        while [ ${curTPR} -le ${tprMax} ];
        do
          # Make sure we are in a suitable range
          numPEsPerNode=$(( ${curPPN} * ${curTPR} ))
          if [ ${minPEcountPerNode} -le ${numPEsPerNode} ] && [ ${maxPEcountPerNode} -ge ${numPEsPerNode} ];
          then
            FullScriptName=${testName}/script_${fileID}id_${roundID}round_${curLaunchID}launchID_${curNumNodes}nodes_${curPPN}ppn_${curTPR}tpr.${BatchFileExtension}
            chmod +x ${FullScriptName}
	    ${Batch} ${FullScriptName}
          fi
          curTPR=$(( ${curTPR} * ${tprScaleFactor} ))
        done
        curPPN=$(( ${curPPN} * ${ppnScaleFactor} ))
      done
      curNumNodes=$(( ${curNumNodes} * ${nodeScaleFactor} ))
      listIndex=$(( ${listIndex} + 1 ))
    done
    curLaunchID=$(( ${curLaunchID} + 1 ))
  done
fi
