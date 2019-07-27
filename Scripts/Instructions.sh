### Below: variables set by user


# ********************************************************************************************************************************
# CritterPath - specify full path to Critter repository
#             - do not include a '/' after 'critter'
CritterPath="${HOME}/src/ExternalLibraries/critter"

# ********************************************************************************************************************************
# MachinePath - machine name corresponding to script file in Machines/ directory
MachinePath="porter"

# ********************************************************************************************************************************
# LibraryPaths - list of library names to be built automatically
#              - append more if necessary
LibraryPaths=("camfs" "candmc")

# ********************************************************************************************************************************
# Algorithm tags (5 predefined)
Algorithms=()

# Fill a list of methods to be tested against, then add to Algorithms list
Test1=('camfs/cqr2' 'candmc/bsqr')

Algorithms+=(${Test1})


# ********************************************************************************************************************************
# fileID - base name of the directory inside which all data/scripts will be stored
#        - will appear inside the SCRATCH directory
fileID=benchQR1

# ********************************************************************************************************************************
# roundID - set to '1' unless performing piecewise testing (launching same job separately) to enhance performance reproducibility
roundID=1

# ********************************************************************************************************************************
# minNumNodes - minimum number of nodes needed for any one test
minNumNodes=1

# ********************************************************************************************************************************
# maxNumNodes - maximum number of nodes needed for any one test
maxNumNodes=1

# ********************************************************************************************************************************
# nodeScaleFactor - scaling factor to apply to the number of nodes
nodeScaleFactor=2

# ********************************************************************************************************************************
# ppnScaleFactor - scaling factor to apply to the number of MPI processes per node (ppn)
ppnScaleFactor=8

# ********************************************************************************************************************************
# tprScaleFactor - scaling factor to apply to the number of threads per MPI rank (tpr)
tprScaleFactor=2

# ********************************************************************************************************************************
# NumLaunchesPerBinary - set to '1' unless performing testing to enhance performance reproducibility (by launching same jobs multiple times)
#                      - different from 'roundID' because the former did not launch at the same time, but waited via a separate launch of bench.sh
NumLaunchesPerBinary=1

# ********************************************************************************************************************************
# numTests - number of scaling studies
#          - for example, weak scaling and strong scaling, even across the same variants, constitute separate tests
numTests=${#Algorithms[@]}

# ********************************************************************************************************************************
# numHours,numMinutes,numSeconds - specify the 
numHours=01
numMinutes=00
numSeconds=00
# ********************************************************************************************************************************
# email - specify email address that you'd like job updates to appear
email="hutter2@illinois.edu"

# ********************************************************************************************************************************
# dataType - float[0], double[1], complex<float>[2], complex<double>[3]
#          - only relevant if test file uses an environment variable to specify this
dataType=1

# ********************************************************************************************************************************
# intType - int[0], int64_t[1]
#         - only relevant if test file uses an environment variable to specify this
intType=1

# ********************************************************************************************************************************
# analyzeDecision1 - profile using Critter
analyzeDecision1=1

# ********************************************************************************************************************************
# analyzeDecision2 - profile using TAU
#                  - not currently supported
analyzeDecision2=0

# ********************************************************************************************************************************
# mpiType - specify 'mpi' (unless on Porter, then 'ampi' is available)
mpiType=mpi

# ********************************************************************************************************************************
# minPEcountPerNode/maxPEcountPerNodempiType - specify the min and max number of processing elements (processes per node x threads per process)
minPEcountPerNode=""
maxPEcountPerNode=""
if [ "$(hostname |grep "porter")" != "" ];
then
  machineName=PORTER
  if [ "${mpiType}" == "mpi" ];
  then
    minPEcountPerNode=1
    maxPEcountPerNode=16
  elif [ "${mpiType}" == "ampi" ];
  then
    minPEcountPerNode=1
    maxPEcountPerNode=512
  fi
elif [ "$(hostname |grep "mira")" != "" ] || [ "$(hostname |grep "cetus")" != "" ];
then
  minPEcountPerNode=16		# Note: this will need to be changed before launching Critter runs
  maxPEcountPerNode=32
elif [ "$(hostname |grep "theta")" != "" ];
then
  minPEcountPerNode=64		# Note: this will need to be changed before launching Critter runs
  maxPEcountPerNode=128
elif [ "$(hostname |grep "stampede2")" != "" ];
then
  minPEcountPerNode=64		# Note: this will need to be changed before launching Critter runs
  maxPEcountPerNode=128
elif [ "$(hostname |grep "h2o")" != "" ];
then
  minPEcountPerNode=16
  maxPEcountPerNode=32
fi

# ********************************************************************************************************************************
# ppnMinList,ppnMaxList,tprMinList,tprMaxList - Lists with an element for each node count
#                                             - specify the min and max ppn,tpr at each node count
#                                             - 'ppn' stands for 'MPI processes per node'
#                                             - 'tpr' stands for 'threads-per-MPI-rank'
ppnMinList=(1)
ppnMaxList=(64)
tprMinList=(1)
tprMaxList=(1)

# ********************************************************************************************************************************
# SubmitToQueue - '1' to submit jobs to queue, '0' to not submit to queue
SubmitToQueue=0
