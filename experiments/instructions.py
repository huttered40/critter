### Below: variables set by user

import os
from operator import (__mul__,__add__)
from critter import bench

from algorithm import (algorithm)
from machines import (porter)
from libraries import (camfs)

CritterPath=os.environ["HOME"]+"/hutter2/ExternalLibraries/critter"
MachineType=porter
LibraryTypeList=[camfs]
fileID="benchQR1"
roundID=1
minNumNodes=1
maxNumNodes=1
nodeScaleFactor=2
ppnScaleFactor=8
tprScaleFactor=2
NumLaunchesPerBinary=1
numTests=1
numHours="01"
numMinutes="00"
numSeconds="00"
email="hutter2@illinois.edu"
dataType=1
intType=1
analyzeDecision1=1
analyzeDecision2=0
mpiType="mpi"
if (os.system("hostname |grep \"porter\"") != ""):
    if (mpiType == "mpi"):
        minPEcountPerNode=1
        maxPEcountPerNode=16
    elif (mpiType == "ampi"):
        minPEcountPerNode=1
        maxPEcountPerNode=512
elif (os.system("hostname |grep \"stampede2\"") != ""):
  minPEcountPerNode=64          # Note: this will need to be changed before launching Critter runs
  maxPEcountPerNode=128
elif (os.system("hostname |grep \"h2o\"") != ""):
  minPEcountPerNode=16
  maxPEcountPerNode=32
ppnMinList=[1]
ppnMaxList=[8]
tprMinList=[1]
tprMaxList=[1]
SubmitToQueue=0
Algorithm1 = algorithm("camfs/cacqr2",\
                       [1024,128,0,0,0,1,1,3],\
		       [1024,128,0,0,0,1,1,3],\
		       [1,1,1,1,1,1,1,1],\
		       [__mul__,__mul__,__mul__,__mul__,__mul__,__mul__,__mul__,__mul__],\
		       [1,1,1,1,1,8,1,1],\
		       [__mul__,__mul__,__mul__,__mul__,__mul__,__mul__,__mul__,__mul__])
Test1=[Algorithm1]
AlgorithmList=[Test1]

Launcher = bench(CritterPath,MachineType,LibraryTypeList,fileID,roundID,minNumNodes,maxNumNodes,nodeScaleFactor,ppnScaleFactor,tprScaleFactor,NumLaunchesPerBinary,\
                 numTests,numHours,numMinutes,numSeconds,email,dataType,intType,analyzeDecision1,analyzeDecision2,mpiType,minPEcountPerNode,maxPEcountPerNode,\
		 ppnMinList,ppnMaxList,tprMinList,tprMaxList,SubmitToQueue,AlgorithmList)
