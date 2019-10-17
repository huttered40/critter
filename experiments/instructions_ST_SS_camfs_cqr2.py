### Below: variables set by user

import os
import math
from operator import (__mul__,__add__,__div__,__sub__)
from critter import (bench,algorithm)

from machines import (stampede2)
from libraries import (camfs)

CritterPath=os.environ["HOME"]+"/critter/"
MachineType=stampede2
LibraryTypeList=[camfs]
CritterVizInfo=[1,1,1048576,4,16,4]
fileID="commcost_vs_ppn"
roundID=1
NumLaunchesPerBinary=1
numTests=1
numHours="00"
numMinutes="30"
numSeconds="00"
email="hutter2@illinois.edu"
minPEcountPerNode=64          # Note: this will need to be changed before launching Critter runs
maxPEcountPerNode=64
nodeMinList=[8]
nodeMaxList=[64]
ppnMinList=[1]
ppnMaxList=[64]
tprMinList=[1]
tprMaxList=[64]
nodeScaleFactorList=[2]
ppnScaleFactorList=[8]
tprScaleFactorList=[8]
nodeScaleOperatorList=[__mul__]
ppnScaleOperatorList=[__mul__]
tprScaleOperatorList=[__mul__]
Algorithm1 = algorithm("camfs_cacqr2",\
                       [16384,2048,1,0,3],\
		       [16384,2048,8,0,3],\
		       [1,1,2,1,1],\
		       [__mul__,__mul__,__mul__,__mul__,__mul__],\
                       lambda x: 0,\
                       lambda InputList,HardwareList: ((((HardwareList[0]*HardwareList[1])/(InputList[2]**2))>=InputList[2]) and (InputList[3] <= int(math.log(InputList[2])))),\
		       [[1,1,1,1,1]],\
		       [[__mul__,__mul__,__mul__,__mul__,__mul__]])
File1 = ["Performance","Residual","Deviation from Orthogonality"]
Test1=[[Algorithm1],"Strong Scaling: 16384x2048 matrix",File1]
TestList=[Test1]

Launcher = bench(CritterPath,MachineType,LibraryTypeList,CritterVizInfo,fileID,roundID,NumLaunchesPerBinary,\
                 numTests,numHours,numMinutes,numSeconds,email,minPEcountPerNode,maxPEcountPerNode,\
		 nodeMinList,nodeMaxList,ppnMinList,ppnMaxList,tprMinList,tprMaxList,nodeScaleFactorList,ppnScaleFactorList,tprScaleFactorList,\
                 nodeScaleOperatorList,ppnScaleOperatorList,tprScaleOperatorList,TestList)
Launcher.generate()
#Launcher.launch()
