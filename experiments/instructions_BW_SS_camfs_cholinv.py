### Below: variables set by user

import os
import math
from operator import (__mul__,__add__,__div__,__sub__)
from critter import (bench,algorithm)

from machines import (bluewaters)
from libraries import (camfs)

CritterPath=os.environ["HOME"]+"/critter/"
MachineType=bluewaters
LibraryTypeList=[camfs]
CritterVizInfo=[1,1,1048576,4,16,4]
fileID="benchCF"
roundID=1
NumLaunchesPerBinary=1
numTests=1
numHours="01"
numMinutes="00"
numSeconds="00"
email="hutter2@illinois.edu"
minPEcountPerNode=16
maxPEcountPerNode=32
nodeMinList=[2]
nodeMaxList=[32]
ppnMinList=[16]
ppnMaxList=[32]
tprMinList=[1]
tprMaxList=[1]
nodeScaleFactorList=[2]
ppnScaleFactorList=[2]
tprScaleFactorList=[2]
nodeScaleOperatorList=[__mul__]
ppnScaleOperatorList=[__mul__]
tprScaleOperatorList=[__mul__]
Algorithm1 = algorithm("camfs_cholinv",\
                       [1024,1,0,3],\
		       [1024,1,0,3],\
		       [1,2,1,1],\
		       [__mul__,__mul__,__mul__,__mul__],\
                       lambda x: 0,\
                       lambda InputList,HardwareList: ((int(round((HardwareList[0]*HardwareList[1])**(1./3.)))**3 == (HardwareList[0]*HardwareList[1])) and (InputList[2] <= round((HardwareList[0]*HardwareList[1])**(1./3.)))),\
		       [[1,1,1,1]],\
		       [[__mul__,__mul__,__mul__,__mul__]])
File1 = ["Performance","Residual"]
Test1=[[Algorithm1],"Strong Scaling",File1]
TestList=[Test1]

Launcher = bench(CritterPath,MachineType,LibraryTypeList,CritterVizInfo,fileID,roundID,NumLaunchesPerBinary,\
                 numTests,numHours,numMinutes,numSeconds,email,minPEcountPerNode,maxPEcountPerNode,\
		 nodeMinList,nodeMaxList,ppnMinList,ppnMaxList,tprMinList,tprMaxList,nodeScaleFactorList,ppnScaleFactorList,tprScaleFactorList,\
                 nodeScaleOperatorList,ppnScaleOperatorList,tprScaleOperatorList,TestList)
Launcher.generate()
Launcher.launch()
