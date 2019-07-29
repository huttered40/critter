class stampede2(object):
    """
    """
    BatchFileExtension="sh"
    Batch="sbatch --mail-user=hutter2@illinois.edu --mailtype=all"
    AllocationName=""
    machineName=STAMPEDE2

    @staticmethod
    def set():
        # Note: SCRATCH environment variable does not need to be set
        pass

    @staticmethod
    def script(scriptFile,testName,curNumNodes,curPPN,curTPR,numPEsPerNode,numHours,numMinutes,numSeconds):
        scriptFile.write("#!/bin/bash\n")
        scriptFile.write("#SBATCH -J %s_%dnodes_%dppn_%dtpr\n" %(testName,curNumNodes,curPPN,curTPR))
        scriptFile.write("#SBATCH -o %s_%dnodes_%dppn_%dtpr.o%j\n" %(testName,curNumNodes,curPPN,curTPR))
        scriptFile.write("#SBATCH -e %s_%dnodes_%dppn_%dtpr.e%j\n" %(testName,curNumNodes,curPPN,curTPR))
        if (curNumNodes <= 256):
            scriptFile.write("#SBATCH -p normal\n")
        else:
            scriptFile.write("#SBATCH -p large\n")
        scriptFile.write("#SBATCH -N %d\n" %(curNumNodes))
        scriptFile.write("#SBATCH -n %d\n" %(curNumNodes*curPPN))
        scriptFile.write("#SBATCH -t %s:%s:%s\n" %(numHours,numMinutes,numSeconds))
        scriptFile.write("export MKL_NUM_THREADS=%d\n" %(curTPR))

    @staticmethod
    def writeTest(scriptFile,numProcesses,ppn,tpr,AlgInputString):
        """
	"""
        Str1="ibrun "
        ScriptFile.write(Str1+AlgInputString)

    @staticmethod
    def queue(Script):
        call("cd %s; chmod +x %s"%(os.environ["SCRATCH"],Script),shell=True)
        call("cd %s; %s %s"%(os.environ["SCRATCH"],Batch,Script),shell=True)
