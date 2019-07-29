### Below: variables set by user

import os
from operator import (__mul__,__add__)
from critter import (bench,algorithm)

from machines import (porter)
from libraries import (camfs)

CritterPath=os.environ["HOME"]+"/hutter2/ExternalLibraries/critter"
MachineType=porter
LibraryTypeList=[camfs]
fileID="benchQR1"
roundID=1
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
nodeMinList=[1]
nodeMaxList=[1]
ppnMinList=[[1]]
ppnMaxList=[[8]]
tprMinList=[[1]]
tprMaxList=[[1]]
nodeScaleFactorList=[2]
ppnScaleFactorList=[8]
tprScaleFactorList=[2]
SubmitToQueue=0
Algorithm1 = algorithm("camfs/cacqr2",\
                       [1024,128,0,0,0,1,1,3],\
		       [1024,128,0,0,0,1,1,3],\
		       [1,1,1,1,1,1,1,1],\
		       [__mul__,__mul__,__mul__,__mul__,__mul__,__mul__,__mul__,__mul__],\
                       .. need lambda to map scaleIndex (0,1,2,3,..) to index into the lists below ..
		       [[1,1,1,1,1,8,1,1]],\
		       [[__mul__,__mul__,__mul__,__mul__,__mul__,__mul__,__mul__,__mul__]])
Test1=[Algorithm1]
AlgorithmList=[Test1]

Launcher = bench(CritterPath,MachineType,LibraryTypeList,fileID,roundID,NumLaunchesPerBinary,\
                 numTests,numHours,numMinutes,numSeconds,email,dataType,intType,analyzeDecision1,analyzeDecision2,mpiType,minPEcountPerNode,maxPEcountPerNode,\
		 nodeMinList,nodeMaxList,ppnMinList,ppnMaxList,tprMinList,tprMaxList,nodeScaleFactorList,ppnScaleFactorList,tprScaleFactorList,SubmitToQueue,AlgorithmList)
