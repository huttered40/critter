### Below: variables set by user

import os
import math
from operator import (__mul__,__add__,__div__,__sub__)
from critter import (bench,algorithm)

from machines import (bluewaters)
from libraries import (candmc)

CritterPath=os.environ["HOME"]+"/critter/"
MachineType=bluewaters
LibraryTypeList=[candmc]
fileID="benchQR"
roundID=1
NumLaunchesPerBinary=1
numTests=1
numHours="00"
numMinutes="30"
numSeconds="00"
email="hutter2@illinois.edu"
minPEcountPerNode=16
maxPEcountPerNode=32
nodeMinList=[4]
nodeMaxList=[32]
ppnMinList=[16]
ppnMaxList=[32]
tprMinList=[1]
tprMaxList=[4]
nodeScaleFactorList=[2]
ppnScaleFactorList=[2]
tprScaleFactorList=[2]
nodeScaleOperatorList=[__mul__]
ppnScaleOperatorList=[__mul__]
tprScaleOperatorList=[__mul__]
Algorithm1 = algorithm("candmc_bench_scala_qr",\
                       [32768,512,8,1,5,0,0,0],\
		       [32768,512,16,4,5,0,0,0],\
		       [1,1,2,2,1,1,1,1],\
		       [__mul__,__mul__,__mul__,__mul__,__mul__,__mul__,__mul__,__mul__],\
                       lambda x: 0,\
                       lambda InputList,HardwareList: (((InputList[1]/InputList[3])>InputList[2]) and ((InputList[0]/(HardwareList[0]*HardwareList[1]/InputList[3])) >= (InputList[1]/InputList[3]))),\
		       [[1,1,1,1,1,1,1,1]],\
		       [[__mul__,__mul__,__mul__,__mul__,__mul__,__mul__,__mul__,__mul__]])
File1 = ["Performance/Node"]
Test1=[[Algorithm1],"Strong Scaling: 8192x512 matrix",File1]
TestList=[Test1]

Launcher = bench(CritterPath,MachineType,LibraryTypeList,fileID,roundID,NumLaunchesPerBinary,\
                 numTests,numHours,numMinutes,numSeconds,email,minPEcountPerNode,maxPEcountPerNode,\
		 nodeMinList,nodeMaxList,ppnMinList,ppnMaxList,tprMinList,tprMaxList,nodeScaleFactorList,ppnScaleFactorList,tprScaleFactorList,\
                 nodeScaleOperatorList,ppnScaleOperatorList,tprScaleOperatorList,TestList)
Launcher.generate()
Launcher.launch()
