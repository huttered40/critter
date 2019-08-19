import os
from subprocess import call

class porter(object):
    """
    """
    BatchFileExtension=""
    Batch=""
    AllocationName=""
    MachineName="PORTER"
    PeakNetworkInjectionRate=None
    PeakNodePerformance=None

    @staticmethod
    def set():
        """
	"""
        os.environ["SCRATCH"] = os.environ["HOME"]+"/hutter2/ExternalLibraries/critter/Tests"
        os.environ["MPITYPE"] = "MPI_TYPE"
        #os.environ["MPITYPE"] = "AMPI_TYPE"

    @staticmethod
    def script(scriptFile,testName,curNumNodes,curPPN,curTPR,numPEsPerNode,numHours,numMinutes,numSeconds):
        """
	"""
	pass

    @staticmethod
    def write_test(ScriptFile,numProcesses,ppn,tpr,AlgInputString):
        """
	"""
        Str1="mpiexec -n %d " %(numProcesses)
        print("\n      "+Str1+AlgInputString)
        call(Str1+AlgInputString,shell=True)
        #${BINARYPATH}charmrun +p1 +vp${numProcesses} ${@:5:$#}

    @staticmethod
    def queue(Script):
        pass

    @staticmethod
    def IsAccelerated():
        return 0
