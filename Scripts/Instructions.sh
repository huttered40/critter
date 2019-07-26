### Below: variables set by user

# ********************************************************************************************************************************
# LibraryPaths - list of paths to libraries to be built automatically
#              - append more if necessary
LibraryPaths=()
if [ "$(hostname |grep "porter")" != "" ];
then
  camfsDir=~/hutter2/CAMFS
  candmcDir=~/hutter2/ExternalLibraries/CANDMC
elif [ "$(hostname |grep "mira")" != "" ] || [ "$(hostname |grep "cetus")" != "" ];
then
  camfsDir=~/scratch/CAMFS
  candmcDir=~/scratch/CANDMC
elif [ "$(hostname |grep "theta")" != "" ];
then
  camfsDir=~/scratch/CAMFS
  candmcDir=~/scratch/CANDMC
elif [ "$(hostname |grep "stampede2")" != "" ];
then
  camfsDir=~/CAMFS
  candmcDir=~/CANDMC
elif [ "$(hostname |grep "h2o")" != "" ];
then
  camfsDir=~/CAMFS
  candmcDir=~/CANDMC
fi
LibraryPaths+=(${camfsDir})
LibraryPaths+=(${candmcDir})

# ********************************************************************************************************************************
# BinaryPath - where binaries (of all libraries) are stored before transfer to SCRATCH
BinaryPath=${LibraryPath}/src/bin

# ********************************************************************************************************************************
# Algorithm tags (5 predefined)
Algorithms=()

Test1=()

# Fill a list of methods to be tested against, then add to Algorithms list
Test1+=('camfs/cqr2')
Test1+=('candmc/bsqr')

Algorithms+=(${Test1})


# ********************************************************************************************************************************
fileID=benchQR1
roundID=1
minNumNodes=1
maxNumNodes=1
nodeScaleFactor=2
ppnScaleFactor=8
tprScaleFactor=2
NumLaunchesPerBinary=1
numTests=${#Algorithms[@]}

numHours=01
numMinutes=00
numSeconds=00
MyEmail="hutter2@illinois.edu"

# ********************************************************************************************************************************
# dataType - float[0], double[1], complex<float>[2], complex<double>[3]
dataType=1

# ********************************************************************************************************************************
# intType - int[0], int64_t[1]
intType=1

# ********************************************************************************************************************************
# analyzeDecision1 - profile using Critter
analyzeDecision1=1

# ********************************************************************************************************************************
# analyzeDecision2 - profile using TAU
analyzeDecision2=0
