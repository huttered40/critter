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
  # Check if we want GPU acceleration (XK7 vs. XE6 nodes)
  #read -p "XE6 (Y) or XK7 (N) node: " nodeType
  #if [ "${nodeType}" == "Y" ];
  #then
  #echo "#PBS -l nodes=${curNumNodes}:ppn=${curPPN}:xe" >> ${scriptName}
  #elif [ "${nodeType}" == "N" ];
  #then
  #  echo "#PBS -l nodes=${curNumNodes}:ppn=${curPPN}:xk" >> ${scriptName}
  #fi
  echo "#PBS -l nodes=${curNumNodes}:ppn=${numPEsPerNode}:xe" >> ${scriptName}
  echo "#PBS -l walltime=${numHours}:${numMinutes}:${numSeconds}" >> ${scriptName}
  echo "#PBS -N camfs" >> ${scriptName}
  echo "#PBS -e ${testName}_${curNumNodes}nodes_${curPPN}ppn_${curTPR}tpr.err" >> ${scriptName}
  echo "#PBS -o ${testName}_${curNumNodes}nodes_${curPPN}ppn_${curTPR}tpr.out" >> ${scriptName}
  echo "##PBS -m Ed" >> ${scriptName}
  echo "#PBS -M ${email}" >> ${scriptName}
  echo "#PBS -A bahv" >> ${scriptName}
  echo "#PBS -W umask=0027" >> ${scriptName}
  #echo "cd ${PBS_O_WORKDIR}" >> ${scriptName}
  echo "#module load craype-hugepages2M  perftools" >> ${scriptName}
  echo "#export APRUN_XFER_LIMITS=1  # to transfer shell limits to the executable" >> ${scriptName}
  echo "export OMP_NUM_THREADS=${curTPR}" >> ${scriptName}
  #if [ "${nodeType}" == "N" ];
  #then
  #  export CRAY_CUDA_MPS=1
  #  export MPICH_RDMA_ENABLED_CUDA=1
  #fi
}
