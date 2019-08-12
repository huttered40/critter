import os
from subprocess import call

class camfs(object):
    """
    """
    @staticmethod
    def build(CritterPath,testName):
        if (os.system("hostname |grep \"porter\"") != 256):
            camfsDir="~/hutter2/CAMFS"
        elif (os.system("hostname |grep \"stampede2\"") != 256):
            camfsDir="~/CAMFS"
        elif (os.system("hostname |grep \"h2o\"") != 256):
            camfsDir="~/CAMFS"

        #call("make -C%s/src clean"%(camfsDir),shell=True)
        #call("make -C%s/src all"%(camfsDir),shell=True)
        #call("cd %s/src/bin; for j in *; do mv -- \"$j\" \"camfs_$j\"; done; cd -;"%(camfsDir),shell=True)
        #call("mv %s/src/bin/* %s/Tests/%s/bin/"%(camfsDir,CritterPath,testName),shell=True)
