import os
from subprocess import call

class porter(object):
    """
    """
    BatchFileExtension=""
    Batch=""
    machineName="PORTER"

    @staticmethod
    def set():
        """
	"""
        os.environ["SCRATCH"] = os.environ["HOME"]+"/hutter2/ExternalLibraries/critter/Tests"

    @staticmethod
    def script(scriptFile,testName,curNumNodes,curPPN,curTPR,numPEsPerNode,numHours,numMinutes,numSeconds):
        """
	"""
	pass

    @staticmethod
    def writeTest(numProcesses,ppn,tpr,AlgInputString):
        """
	"""
	Str1="mpiexec -n %d " %(numProcesses)
        #print(Str1+AlgInputString)
        call(Str1+AlgInputString,shell=True)
        #${BINARYPATH}charmrun +p1 +vp${numProcesses} ${@:5:$#}

    @staticmethod
    def queue(Script):
        pass

    @staticmethod
    def IsAccelerated():
        return 0
