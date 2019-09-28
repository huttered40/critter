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
fileID="benchQRWS"
roundID=1
NumLaunchesPerBinary=1
numTests=1
numHours="01"
numMinutes="00"
numSeconds="00"
email="hutter2@illinois.edu"
minPEcountPerNode=64          # Note: this will need to be changed before launching Critter runs
maxPEcountPerNode=64
nodeMinList=[8]
nodeMaxList=[128]
ppnMinList=[32]
ppnMaxList=[64]
tprMinList=[1]
tprMaxList=[2]
nodeScaleFactorList=[2]
ppnScaleFactorList=[2]
tprScaleFactorList=[2]
nodeScaleOperatorList=[__mul__]
ppnScaleOperatorList=[__mul__]
tprScaleOperatorList=[__mul__]
Algorithm1 = algorithm("camfs_cacqr2",\
                       [8192,1024,2,0,5],\
		       [8192,1024,4,2,5],\
		       [1,1,2,1,1],\
		       [__mul__,__mul__,__mul__,__add__,__mul__],\
                       lambda x: (1 if (x%4==0) else 0),\
                       lambda InputList,HardwareList: ((((HardwareList[0]*HardwareList[1])/(InputList[2]**2))>=InputList[2]) and (InputList[3] <= int(math.log(InputList[2],2)))),\
		       [[2,1,1,1,1],[2,2,2,1,1]],\
		       [[__mul__,__mul__,__mul__,__mul__,__mul__],[__div__,__mul__,__mul__,__mul__,__mul__]])
File1 = ["Performance","Residual","Deviation from Orthogonality"]
Test1=[[Algorithm1],"Weak Scaling: 8192x512 initial matrix",File1]
TestList=[Test1]

Launcher = bench(CritterPath,MachineType,LibraryTypeList,fileID,roundID,NumLaunchesPerBinary,\
                 numTests,numHours,numMinutes,numSeconds,email,minPEcountPerNode,maxPEcountPerNode,\
		 nodeMinList,nodeMaxList,ppnMinList,ppnMaxList,tprMinList,tprMaxList,nodeScaleFactorList,ppnScaleFactorList,tprScaleFactorList,\
                 nodeScaleOperatorList,ppnScaleOperatorList,tprScaleOperatorList,TestList)
#Launcher.build()
Launcher.generate()
#Launcher.launch()
