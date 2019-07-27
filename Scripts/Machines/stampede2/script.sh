launch () {
  scriptName=${1}
  testName=${2}
  curNumNodes=${3}
  curPPN=${4}
  curTPR={5}
  numPEsPerNode={6}
  numHours={7}
  numMinutes={8}
  numSeconds={9}
  echo "#!/bin/bash" > ${scriptName}
  echo "#SBATCH -J ${testName}_${curNumNodes}nodes_${curPPN}ppn_${curTPR}tpr" >> ${scriptName}
  echo "#SBATCH -o ${testName}_${curNumNodes}nodes_${curPPN}ppn_${curTPR}tpr.o%j" >> ${scriptName}
  echo "#SBATCH -e ${testName}_${curNumNodes}nodes_${curPPN}ppn_${curTPR}tpr.e%j" >> ${scriptName}
  if [ ${curNumNodes} -le 256 ];
  then
    echo "#SBATCH -p normal" >> ${scriptName}
  else
    echo "#SBATCH -p large" >> ${scriptName}
  fi
  echo "#SBATCH -N ${curNumNodes}" >> ${scriptName}
  echo "#SBATCH -n $((${curNumNodes} * ${curPPN}))" >> ${scriptName}
  echo "#SBATCH -t ${numHours}:${numMinutes}:${numSeconds}" >> ${scriptName}
  echo "export MKL_NUM_THREADS=${curTPR}" >> ${scriptName}
}
