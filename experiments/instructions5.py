### Below: variables set by user

import os
import math
from operator import (__mul__,__add__,__div__,__sub__)
from critter import (bench,algorithm)

from machines import (porter)
from libraries import (candmc)

Debug=1
CritterPath=os.environ["HOME"]+"/hutter2/ExternalLibraries/critter"
MachineType=porter
LibraryTypeList=[candmc]
fileID="benchQR1"
roundID=1
NumLaunchesPerBinary=1
numTests=1
numHours="01"
numMinutes="00"
numSeconds="00"
email="hutter2@illinois.edu"
mpiType="mpi"
if (os.system("hostname |grep \"porter\"") != 256):
    if (mpiType == "mpi"):
        minPEcountPerNode=64
        maxPEcountPerNode=128
    elif (mpiType == "ampi"):
        minPEcountPerNode=1
        maxPEcountPerNode=512
elif (os.system("hostname |grep \"stampede2\"") != 256):
    minPEcountPerNode=64          # Note: this will need to be changed before launching Critter runs
    maxPEcountPerNode=128
elif (os.system("hostname |grep \"h2o\"") != 256):
    minPEcountPerNode=16
    maxPEcountPerNode=32
nodeMinList=[8]
nodeMaxList=[4096]
ppnMinList=[[8,8,8,8,8,8,8,8,8,8]]
ppnMaxList=[[64,64,64,64,64,64,64,64,64,64]]
tprMinList=[[1,1,1,1,1,1,1,1,1,1]]
tprMaxList=[[2,2,2,2,2,2,2,2,2,2]]
nodeScaleFactorList=[2]
ppnScaleFactorList=[2]
tprScaleFactorList=[2]
nodeScaleOperatorList=[__mul__]
ppnScaleOperatorList=[__mul__]
tprScaleOperatorList=[__mul__]
Algorithm1 = algorithm("candmc_bsqr",\
                       [8388608,512,16,3,1,0,0],\
		       [8388608,512,64,3,2,0,0],\
		       [1,1,2,1,2,0,0],\
		       [__mul__,__mul__,__mul__,__add__,__mul__,__add__,__add__],\
                       lambda x: (1 if (x%4==0) else 0),\
                       lambda InputList,HardwareList: True,\
		       [[2,1,1,1,1,1,1],[2,2,1,1,2,1,1]],\
		       [[__mul__,__mul__,__mul__,__mul__,__mul__,__mul__,__mul__],[__div__,__mul__,__mul__,__mul__,__mul__,__mul__,__mul__]])
File1 = [["perf",("Performance","Residual","Deviation from Orthogonality")]]
Test1=[[Algorithm1],"Weak Scaling",File1]
TestList=[Test1]

Launcher = bench(CritterPath,MachineType,LibraryTypeList,fileID,roundID,NumLaunchesPerBinary,\
                 numTests,numHours,numMinutes,numSeconds,email,minPEcountPerNode,maxPEcountPerNode,\
		 nodeMinList,nodeMaxList,ppnMinList,ppnMaxList,tprMinList,tprMaxList,nodeScaleFactorList,ppnScaleFactorList,tprScaleFactorList,\
                 nodeScaleOperatorList,ppnScaleOperatorList,tprScaleOperatorList,TestList)
Launcher.build()
#Launcher.generate()
#Launcher.launch()
