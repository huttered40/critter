class candmc(object):
    """
    """
    @staticmethod
    def build():
        if [ "$(hostname |grep "porter")" != "" ];
        then
            candmcDir=~/hutter2/ExternalLibraries/CANDMC
        elif [ "$(hostname |grep "stampede2")" != "" ];
        then
            candmcDir=~/CANDMC
        elif [ "$(hostname |grep "h2o")" != "" ];
        then
            candmcDir=~/CANDMC
        fi

        cd ${candmcDir}
        make clean
        rm config.mk
        ./configure
        make bench_scala_qr
        cd -
        mv ${candmcDir}/bin/benchmarks/bench_scala_qr ${candmcDir}/bin/benchmarks/candmc_rsqr_${machineName}_${PROFTYPE}
        mv ${candmcDir}/bin/benchmarks/candmc_rsqr_${machineName}_${PROFTYPE} ${CritterPath}/Tests/${testName}/bin/
