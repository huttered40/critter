from subprocess import call

class porter(object):
    """
    """
    BatchFileExtension=""
    Batch=""
    machineName=PORTER

    @staticmethod
    def set():
        """
	"""
        os.environ["SCRATCH"] = os.environ["HOME"]+"/hutter2/ExternalLibraries/critter/Tests"

    @staticmethod
    def script():
        """
	"""
	pass

    @staticmethod
    def writeTest(scriptFile,numProcesses,ppn,tpr,AlgInputString):
        """
	"""
	Str1="mpiexec -n %d " %(numProcesses)
        call(Str1+AlgInputString,,shell=True))
        #${BINARYPATH}charmrun +p1 +vp${numProcesses} ${@:5:$#}

    @staticmethod
    def queue(Script):
        pass
        call("chmod +x ${Script}"%(Script),shell=True)
        call("${Batch} ${Script}"%(Batch,Script),shell=True)
