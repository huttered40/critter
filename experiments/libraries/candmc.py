import os
from subprocess import call

class candmc(object):
    """
    """
    @staticmethod
    def build(CritterPath,testName):
        if (os.system("hostname |grep \"porter\"") != 256):
            candmcDir="~/hutter2/ExternalLibraries/candmc"
        elif (os.system("hostname |grep \"stampede2\"") != 256):
            candmcDir="~/candmc"
        elif (os.system("hostname |grep \"h2o\"") != 256):
            candmcDir="~/candmc"

        call("cd %s; make clean; rm config.mk; ./configure; make bench; cd -"%(candmcDir),shell=True)
        call("cd %s/bin/benchmarks/; for j in *; do mv -- \"$j\" \"candmc_$j\"; done; cd -;"%(candmcDir),shell=True)
        call("mv %s/bin/benchmarks/* %s/Tests/%s/bin/"%(candmcDir,CritterPath,testName),shell=True)
