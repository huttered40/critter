# Choice of compiler for Blue Waters (assumes Cray compiler is loaded by default)
if [ "${machineName}" == "BLUEWATERS" ];
then
  read -p "Do you want the Intel Programming Environment (I) or the GNU Programming Environment (G) (choose G if running on GPU): " bwPrgEnv
  if [ "${bwPrgEnv}" == "I" ];
  then
    if [ "${PE_ENV}" == "GNU" ];
    then
      module swap PrgEnv-gnu PrgEnv-intel
    elif [ "${PE_ENV}" == "CRAY" ];
    then
      module swap PrgEnv-cray PrgEnv-intel
    fi
  elif [ "${bwPrgEnv}" == "G" ];
  then
    if [ "${PE_ENV}" == "INTEL" ];
    then
      module swap PrgEnv-intel PrgEnv-gnu
    elif [ "${PE_ENV}" == "CRAY" ];
    then
      module swap PrgEnv-cray PrgEnv-gnu
    fi
  fi
  if [ "${accelType}" == "n" ];
  then
    module load cblas
  else
    module load cudatoolkit
    # Swap or load anything else? Does the PrgEnv matter with Cuda?
  fi
fi
export SCRATCH=/scratch/sciteam/hutter
fileExtension="pbs"
Batch="qsub"
