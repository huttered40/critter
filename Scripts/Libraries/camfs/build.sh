build_camfs () {
if [ "$(hostname |grep "porter")" != "" ];
then
  camfsDir=~/hutter2/CAMFS
elif [ "$(hostname |grep "stampede2")" != "" ];
then
  camfsDir=~/CAMFS
elif [ "$(hostname |grep "h2o")" != "" ];
then
  camfsDir=~/CAMFS
fi

make -C${camfsDir}/src clean
make -C${camfsDir}/src all
mv ${camfsDir}/src/bin/* ${CritterPath}/Tests/${testName}/bin/
}
