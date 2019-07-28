launch () {
  echo "No-op"
}

writeTest () {
  local numProcesses=${1}
  local ppn=${2}
  local tpr=${3}
  local scriptName=${4}
  if [ "${mpiType}" == "mpi" ];
  then
    mpiexec -n ${numProcesses} ${@:5:$#}
  elif [ "${mpiType}" == "ampi" ];
  then
    ${BINARYPATH}charmrun +p1 +vp${numProcesses} ${@:5:$#}
  fi
}
