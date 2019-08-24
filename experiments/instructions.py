### Below: variables set by user

import os
import math
from operator import (__mul__,__add__,__div__,__sub__)
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
mpiType="mpi"	# Not specified to critter
if (os.system("hostname |grep \"porter\"") != 256):
    if (mpiType == "mpi"):
        minPEcountPerNode=1
        maxPEcountPerNode=16
    elif (mpiType == "ampi"):
        minPEcountPerNode=1
        maxPEcountPerNode=512
elif (os.system("hostname |grep \"stampede2\"") != 256):
    minPEcountPerNode=64          # Note: this will need to be changed before launching Critter runs
    maxPEcountPerNode=128
elif (os.system("hostname |grep \"h2o\"") != 256):
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
nodeScaleOperatorList=[__mul__]
ppnScaleOperatorList=[__mul__]
tprScaleOperatorList=[__mul__]
Algorithm1 = algorithm("camfs_cacqr2",\
                       [512,64,1,0,0,0,3],\
		       [512,64,2,0,0,0,3],\
		       [1,1,2,1,1,1,1],\
		       [__mul__,__mul__,__mul__,__mul__,__mul__,__mul__,__mul__],\
                       lambda x: 0,\
                       lambda InputList,HardwareList: ((((HardwareList[1]*HardwareList[2])/(InputList[2]**2))>=InputList[2]) and (InputList[3] <= int(math.log(InputList[2])))),\
		       [[1,1,1,1,1,1,1]],\
		       [[__mul__,__mul__,__mul__,__mul__,__mul__,__mul__,__mul__]],\
                       [0])
File1 = [["critter",()],["perf",("Performance","Residual","Deviation from Orthogonality")]]
Test1=[[Algorithm1],"Strong Scaling",File1]
TestList=[Test1]

Launcher = bench(CritterPath,MachineType,LibraryTypeList,fileID,roundID,NumLaunchesPerBinary,\
                 numTests,numHours,numMinutes,numSeconds,email,minPEcountPerNode,maxPEcountPerNode,\
		 nodeMinList,nodeMaxList,ppnMinList,ppnMaxList,tprMinList,tprMaxList,nodeScaleFactorList,ppnScaleFactorList,tprScaleFactorList,\
                 nodeScaleOperatorList,ppnScaleOperatorList,tprScaleOperatorList,TestList)
Launcher.build()
Launcher.generate()
Launcher.launch()
