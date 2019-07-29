from subprocess import call

class candmc(object):
    """
    """
    @staticmethod
    def build(CritterPath,testName):
        if (os.system("hostname |grep \"porter\"") != ""):
        then
            candmcDir="~/hutter2/ExternalLibraries/CANDMC"
        elif (os.system("hostname |grep \"stampede2\"") != ""):
        then
            candmcDir="~/CANDMC"
        elif (os.system("hostname |grep \"h2o\"") != ""):
        then
            candmcDir="~/CANDMC"
        fi

        call("cd %s; make clean; rm config.mk; ./configure; make bench_scala_qr; cd -"%(candmcDir),shell=True)
        call("mv %s/bin/benchmarks/bench_scala_qr %s/bin/benchmarks/candmc_rsqr_%s"%(candmcDir,candmcDir,os.environ["PROFTYPE"]),shell=True)
        call("mv %s/bin/benchmarks/candmc_rsqr_%s %s/Tests/%s/bin/"%(candmcDir,os.environ["PROFTYPE"],CritterPath,testName),shell=True)
