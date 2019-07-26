#!/bin/bash

# Warning to set relevant paths before proceeding with launch
read -p "Warning: Have you set the user-defined environment variables inside bench.sh? Yes[1] No[0]" FillEnvironVar

source Instructions.sh

# Create Build Instructions for all libraries
for lib in "${LibraryPaths[@]}"
do
  source Libraries/${lib}/build.sh
done

# Make sure that the bin directory is created, or else compilation won't work
if [ ! -d "${BinaryPath}" ];
then
  mkdir ${BinaryPath}
fi

machineName=""
accelType=""
testAccel_NoAccel=""
if [ "$(hostname |grep "porter")" != "" ];
then
  machineName=PORTER
elif [ "$(hostname |grep "mira")" != "" ] || [ "$(hostname |grep "cetus")" != "" ];
then
  machineName=BGQ
elif [ "$(hostname |grep "theta")" != "" ];
then
  machineName=THETA
elif [ "$(hostname |grep "stampede2")" != "" ];
then
  machineName=STAMPEDE2
elif [ "$(hostname |grep "h2o")" != "" ];
then
  accelType="n"
  testAccel_NoAccel="n"
#  read -p "GPU acceleration via XK7[y] or no[n]: " accelType
#  if [ "${accelType}" == "y" ];
#  then
#    export GPU=GPUACCEL
#    read -p "Do you want to test CA-CQR2 on both GPU accelated machines and the non-accelerated option[y] or no[n]: " testAccel_NoAccel
#  else
#    export GPU=NOGPU
#    testAccel_NoAccel="n"
#  fi
  machineName=BLUEWATERS
fi

if [ "${mpiType}" == "mpi" ];
then
  export MPITYPE=MPI_TYPE
elif [ "${mpiType}" == "ampi" ];
then
  export MPITYPE=AMPI_TYPE
fi

dateStr=$(date +%Y-%m-%d-%H_%M_%S)
ppnMinList=()
ppnMaxList=()
tprMinList=()
tprMaxList=()
ppnMin=""
ppnMax=""
tprMin=""
tprMax=""
read -p "Will ppn/tpr be the same for all tested node counts across all tests? yes[1] or no[0]: " skipHardwareSpec
if [ ${skipHardwareSpec} == 1 ];
then
  read -p "Enter min ppn: " ppnMin
  read -p "Enter max ppn: " ppnMax
  # Only revelant for non-GPU
  read -p "Enter min tpr: " tprMin
  read -p "Enter max tpr: " tprMax
fi
curNumNodes=${minNumNodes}
while [ ${curNumNodes} -le ${maxNumNodes} ];
do
  if [ ${skipHardwareSpec} == 0 ];
  then
    read -p "Enter min ppn for node count ${curNumNodes}: " ppnMin
    read -p "Enter max ppn for node count ${curNumNodes}: " ppnMax
    # Only revelant for non-GPU
    read -p "Enter min tpr for node count ${curNumNodes}: " tprMin
    read -p "Enter max tpr for node count ${curNumNodes}: " tprMax
  fi
  ppnMinList+=(${ppnMin})
  ppnMaxList+=(${ppnMax})
  tprMinList+=(${tprMin})
  tprMaxList+=(${tprMax})

  curNumNodes=$(( ${curNumNodes} * ${nodeScaleFactor} ))   # So far, only use cases for nodeScaleFactor are 2 and 16.
done

fileName=launch${fileID}_${dateStr}_${machineName}_round${roundID}
fileNameToProcess=launch${fileID}_${machineName}	# Name of the corresponding directory in CAMFS_data. Allows for appending multiple runs

if [ ${dataType} == 0 ];
then
  export DATATYPE=FLOAT_TYPE
elif [ ${dataType} == 1 ];
then
  export DATATYPE=DOUBLE_TYPE
elif [ ${dataType} == 2 ];
then
  export DATATYPE=COMPLEX_FLOAT_TYPE
elif [ ${dataType} == 3 ];
then
  export DATATYPE=COMPLEX_DOUBLE_TYPE
fi
if [ ${intType} == 0 ];
then
  export INTTYPE=INT_TYPE
elif [ ${intType} == 1 ];
then
  export INTTYPE=INT64_T_TYPE
fi

###################################################### Library Builds ######################################################

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

# Build each library
for lib in "${LibraryPaths[@]}"
do
  export PROFTYPE=PERFORMANCE
  profType=P
  # export SPECIAL_SCALA_ARG=REF
  build_${lib}

  if [ ${analyzeDecision1} == 1 ];
  then
    profType=${profType}C
    export PROFTYPE=CRITTER
    make -C./.. ${makebinarytag}_${mpiType}
  fi

  if [ ${analyzeDecision2} == 1 ];
  then
    profType=${profType}T
    export PROFTYPE=PROFILE
    make -C./.. ${makebinarytag}_${mpiType}
  fi
done

if [ "${machineName}" == "BGQ" ];
then
  export SCRATCH=/projects/QMCat/huttered
elif [ "${machineName}" == "THETA" ];
then
  export SCRATCH=/projects/QMCat/huttered
  export BINPATH=${SCRATCH}/${fileName}/bin/
elif [ "${machineName}" == "STAMPEDE2" ];
then
  export BINPATH=${SCRATCH}/${fileName}/bin/
elif [ "${machineName}" == "BLUEWATERS" ];
then
  export SCRATCH=/scratch/sciteam/hutter
  export BINPATH=${SCRATCH}/${fileName}/bin/
elif [ "${machineName}" == "PORTER" ];
then
  export SCRATCH=${HOME}/hutter2/Critter_data
  export BINPATH=${BinaryPath}
fi

# collectData.sh will always be a single line, just a necessary intermediate step.
echo "bash $SCRATCH/${fileName}/collectInstructionsStage1.sh | bash PackageDataRemoteStage1.sh" > collectData.sh

cat <<-EOF > $SCRATCH/${fileName}.sh
scriptName=$SCRATCH/${fileName}/script.sh
mkdir $SCRATCH/${fileName}/
mkdir $SCRATCH/${fileName}/DataFiles

# Need to re-build ppn/tpr lists (for each node count) because I cannot access the pre-time list with run-time indices
ppnMinListRunTime=()
ppnMaxListRunTime=()
tprMinListRunTime=()
tprMaxListRunTime=()

ppnMin=""
ppnMax=""
tprMin=""
tprMax=""
read -p "Will ppn/tpr be the same for all tested node counts across all tests? yes[1] or no[0]: " skipHardwareSpec
if [ \${skipHardwareSpec} == 1 ];
then
  read -p "Enter min ppn: " ppnMin
  read -p "Enter max ppn: " ppnMax
  # Only revelant for non-GPU
  read -p "Enter min tpr: " tprMin
  read -p "Enter max tpr: " tprMax
fi
curNumNodes=${minNumNodes}
while [ \${curNumNodes} -le ${maxNumNodes} ];
do
  if [ \${skipHardwareSpec} == 0 ];
  then
    read -p "Enter min ppn for node count \${curNumNodes}: " ppnMin
    read -p "Enter max ppn for node count \${curNumNodes}: " ppnMax
    # Only revelant for non-GPU
    read -p "Enter min tpr for node count \${curNumNodes}: " tprMin
    read -p "Enter max tpr for node count \${curNumNodes}: " tprMax
  fi
  ppnMinListRunTime+=(\${ppnMin})
  ppnMaxListRunTime+=(\${ppnMax})
  tprMinListRunTime+=(\${tprMin})
  tprMaxListRunTime+=(\${tprMax})

  curNumNodes=\$(( \${curNumNodes} * ${nodeScaleFactor} ))   # So far, only use cases for nodeScaleFactor are 2 and 16.
done

curLaunchID=1
while [ \${curLaunchID} -le ${NumLaunchesPerBinary} ];
do
  # Loop over all scripts - log(P) of them
  curNumNodes=${minNumNodes}
  ppnIndex=0
  while [ \${curNumNodes} -le ${maxNumNodes} ];
  do
    minPPN=\${ppnMinListRunTime[\${ppnIndex}]}
    maxPPN=\${ppnMaxListRunTime[\${ppnIndex}]}
    curPPN=\${minPPN}
    while [ \${curPPN} -le \${maxPPN} ];
    do
      minTPR=\${tprMinListRunTime[\${ppnIndex}]}
      maxTPR=\${tprMaxListRunTime[\${ppnIndex}]}
      curTPR=\${minTPR}
      while [ \${curTPR} -le \${maxTPR} ];
      do
        # Make sure we are in a suitable range
        numPEsPerNode=\$(( \${curPPN} * \${curTPR} ))
        if [ ${minPEcountPerNode} -le \${numPEsPerNode} ] && [ ${maxPEcountPerNode} -ge \${numPEsPerNode} ];
        then
	  scriptName=$SCRATCH/${fileName}/script_${fileID}id_${roundID}round_\${curLaunchID}launchID_\${curNumNodes}nodes_\${curPPN}ppn_\${curTPR}tpr
          if [ "${machineName}" == "BGQ" ];
          then
            scriptName=\${scriptName}.sh
            echo "TODO" > \${scriptName}
          fi
	  if [ "${machineName}" == "BLUEWATERS" ];
	  then
            scriptName=\${scriptName}.pbs
	    echo "#!/bin/bash" > \${scriptName}
            # Check if we want GPU acceleration (XK7 vs. XE6 nodes)
            #read -p "XE6 (Y) or XK7 (N) node: " nodeType
            #if [ "${nodeType}" == "Y" ];
            #then
	    #echo "#PBS -l nodes=\${curNumNodes}:ppn=\${curPPN}:xe" >> \${scriptName}
            #elif [ "${nodeType}" == "N" ];
            #then
	    #  echo "#PBS -l nodes=\${curNumNodes}:ppn=\${curPPN}:xk" >> \${scriptName}
            #fi
	    echo "#PBS -l nodes=\${curNumNodes}:ppn=\${numPEsPerNode}:xe" >> \${scriptName}
	    echo "#PBS -l walltime=${numHours}:${numMinutes}:${numSeconds}" >> \${scriptName}
	    echo "#PBS -N camfs" >> \${scriptName}
	    echo "#PBS -e ${fileName}_\${curNumNodes}_\${curPPN}.err" >> \${scriptName}
	    echo "#PBS -o ${fileName}_\${curNumNodes}_\${curPPN}.out" >> \${scriptName}
	    echo "##PBS -m Ed" >> \${scriptName}
	    echo "#PBS -M hutter2@illinois.edu" >> \${scriptName}
	    echo "#PBS -A bahv" >> \${scriptName}
	    echo "#PBS -W umask=0027" >> \${scriptName}
  #          echo "cd \${PBS_O_WORKDIR}" >> \${scriptName}
	    echo "#module load craype-hugepages2M  perftools" >> \${scriptName}
	    echo "#export APRUN_XFER_LIMITS=1  # to transfer shell limits to the executable" >> \${scriptName}
            echo "export OMP_NUM_THREADS=\${curTPR}" >> \${scriptName}
            #if [ "${nodeType}" == "N" ];
            #then
            #  export CRAY_CUDA_MPS=1
            #  export MPICH_RDMA_ENABLED_CUDA=1
            #fi
	  elif [ "${machineName}" == "THETA" ];
	  then
            scriptName=\${scriptName}.sh
	    echo "#!/bin/bash" > \${scriptName}
	    echo "#COBALT -t ${numMinutes}" >> \${scriptName}
	    echo "#COBALT -n \${curNumNodes}" >> \${scriptName}
	    echo "#COBALT --attrs mcdram=cache:numa=quad" >> \${scriptName}
	    echo "#COBALT -A QMCat" >> \${scriptName}
	    echo "export n_nodes=\${curNumNodes}" >> \${scriptName}
	    echo "export n_mpi_ranks_per_node=\${curPPN}" >> \${scriptName}
	    echo "export n_mpi_ranks=\$((\${curNumNodes} * \${curPPN}))" >> \${scriptName}
	    read -p "Enter number of OpenMP threads per rank: " numOMPthreadsPerRank
	    read -p "Enter number of hyperthreads per core: " numHyperThreadsPerCore
	    read -p "Enter number of hyperthreads skipped per rank: " numHyperThreadsSkippedPerRank
	    echo "export n_openmp_threads_per_rank=\${numOMPthreadsPerRank}" >> \${scriptName}
	    echo "export n_hyperthreads_per_core=\${numHyperThreadsPerCore}" >> \${scriptName}
	    echo "export n_hyperthreads_skipped_between_ranks=\${numHyperThreadsSkippedPerRank}" >> \${scriptName}
	  elif [ "${machineName}" == "STAMPEDE2" ];
	  then
            scriptName=\${scriptName}.sh
	    echo "bash script name: \${scriptName}"
	    echo "#!/bin/bash" > \${scriptName}
	    echo "#SBATCH -J myjob_${fileID}id_${roundID}round_\${curNumNodes}nodes_\${curPPN}ppn_\${curTPR}tpr" >> \${scriptName}
	    echo "#SBATCH -o myjob_${fileID}id_${roundID}round_\${curNumNodes}nodes_\${curPPN}ppn_\${curTPR}tpr.o%j" >> \${scriptName}
	    echo "#SBATCH -e myjob_${fileID}id_${roundID}round_\${curNumNodes}nodes_\${curPPN}ppn_\${curTPR}tpr.e%j" >> \${scriptName}
	    if [ \${curNumNodes} -le 256 ];
	    then
	      echo "#SBATCH -p normal" >> \${scriptName}
	    else
	      echo "#SBATCH -p large" >> \${scriptName}
	    fi
	    echo "#SBATCH -N \${curNumNodes}" >> \${scriptName}
	    echo "#SBATCH -n \$((\${curNumNodes} * \${curPPN}))" >> \${scriptName}
	    echo "#SBATCH -t ${numHours}:${numMinutes}:${numSeconds}" >> \${scriptName}
	    echo "export MKL_NUM_THREADS=\${curTPR}" >> \${scriptName}
	  fi
        fi
        curTPR=\$(( \${curTPR} * ${tprScaleFactor} ))
      done
      curPPN=\$(( \${curPPN} * ${ppnScaleFactor} ))
    done
    curNumNodes=\$(( \${curNumNodes} * ${nodeScaleFactor} ))   # So far, only use cases for nodeScaleFactor are 2 and 16.
    ppnIndex=\$(( \${ppnIndex} + 1 ))
  done
  curLaunchID=\$(( \${curLaunchID} + 1 ))
done

# Now I need to use a variable for the command-line prompt, since it will change based on the binary executable,
#   for example, scalapack QR has multiple special inputs that take up comm-line prompts that others dont
#   I want special functions in the inner-loop to handle this

updateCounter () {
  local counter=\${1}
  if [ \${2} -eq 1 ];
  then
    counter=\$((\${counter} + \${3})) 
  elif [ \${2} -eq 2 ];
  then
   counter=\$((\${counter} - \${3})) 
  elif [ \${2} -eq 3 ];
  then
    counter=\$((\${counter} * \${3})) 
  elif [ \${2} -eq 4 ];
  then
    counter=\$((\${counter} / \${3})) 
  fi
  echo "\${counter}"
}

# Note: this function is only used for finding the number of dependencies for each binary run
# Therefore, I can multiply the output by 2 (total and average) and it won't affect anything else
# New note: I am getting rid of the *2, since that will be understood by all other scripts, including MakePlotScript.sh
findCountLength () {
  local curr=\${1}
  local counter=0
  while [ \${curr} -le \${2} ];
  do
    curr=\$(updateCounter \${curr} \${3} \${4})
    counter=\$(( counter+1 ))
  done
  echo "\${counter}"
}

log2 () {
    local x=0
    for (( y=\${1}-1 ; \${y} > 0; y >>= 1 )) ; do
        let x=\${x}+1
    done
    echo \${x}
}


# Writes the beginning of collectInstructionsStage1 and collectInstructionsStage2
WriteHeaderForCollection () {
  echo "echo \"${fileName}\"" > \${1}
  echo "echo \"${fileNameToProcess}\"" >> \${1}
  echo "echo \"${machineName}\"" >> \${1}
  echo "echo \"${profType}\"" >> \${1}
  echo "echo \"${nodeScaleFactor}\"" >> \${1}
  echo "echo \"${numTests}\"" >> \${1}
}


# Non-scalapack
writePlotFileName() {
  # Performance runs will always run, so no reason for an if-statement here
  Prefix1=""
  Prefix2=""
  if [ "\${3}" == "1" ];
  then
    Prefix1="Raw/"
    Prefix2="Stats/"
  fi
  echo "echo \"\${Prefix1}\${1}_perf.txt\"" >> \${2}
  if [ "\${3}" == "1" ];
  then
    echo "echo \"\${Prefix2}\${1}_perf_stats.txt\"" >> \${2}
  fi
  
  echo "echo \"\${Prefix1}\${1}_numerics.txt\"" >> \${2}
  if [ "\${3}" == "1" ];
  then
    echo "echo \"\${Prefix2}\${1}_numerics_stats.txt\"" >> \${2}
  fi

  if [ "${profType}" == "PC" ] || [ "${profType}" == "PCT" ];
  then
    echo "echo \"\${1}_critter.txt\"" >> \${2}
    if [ "\${3}" == "1" ];
    then
      echo "echo \"\${1}_critter_breakdown.txt\"" >> \${2}
    fi
  fi
  if [ "${profType}" == "PT" ] || [ "${profType}" == "PCT" ];
  then
    echo "echo \"\${1}_timer.txt\"" >> \${2}
  fi
}

# Only for bsqr/rsqr -- only necessary for Performance now. Might want to use Critter later, but not Profiler
writePlotFileNameScalapackQR() {
  Prefix1=""
  Prefix2=""
  if [ "\${3}" == "1" ];
  then
    Prefix1="Raw/"
    Prefix2="Stats/"
  fi
  echo "echo \"\${Prefix1}\${1}_NoFormQ.txt\"" >> \${2}
  if [ "\${3}" == "1" ];
  then
    echo "echo \"\${Prefix2}\${1}_NoFormQ_stats.txt\"" >> \${2}
  fi
  echo "echo \"\${Prefix1}\${1}_FormQ.txt\"" >> \${2}
  if [ "\${3}" == "1" ];
  then
    echo "echo \"\${Prefix2}\${1}_FormQ_stats.txt\"" >> \${2}
  fi
}

# Only for bscf -- only necessary for Performance now. Might want to use Critter later, but not Profiler
writePlotFileNameScalapackCholesky() {
  Prefix1=""
  Prefix2=""
  if [ "\${3}" == "1" ];
  then
    Prefix1="Raw/"
    Prefix2="Stats/"
  fi
  echo "echo \"\${Prefix1}\${1}.txt\"" >> \${2}
  if [ "\${3}" == "1" ];
  then
    echo "echo \"\${Prefix2}\${1}_stats.txt\"" >> \${2}
  fi
}


# Functions that write the actual script, depending on machine
launchJobs () {
  local launchID=\${3}
  local numNodes=\${4}
  local ppn=\${5}
  local tpr=\${6}
  local numProcesses=\$((\${numNodes} * \${ppn}))
  local scriptName=script_${fileID}id_${roundID}round_\${launchID}launchID_\${numNodes}nodes_\${ppn}ppn_\${tpr}tpr
  echo "What is scriptName - \${scriptName}"
  if [ "$machineName" == "BGQ" ];
  then
    echo "runjob --np \${numProcesses} -p \${ppn} --block \$COBALT_PARTNAME --verbose=INFO : \${@:7:\$#}" >> $SCRATCH/${fileName}/\${scriptName}.sh
  elif [ "$machineName" == "BLUEWATERS" ];
  then
    echo "aprun -n \${numProcesses} -N \${ppn} -d \${tpr} \${@:7:\$#}" >> $SCRATCH/${fileName}/\${scriptName}.pbs
  elif [ "$machineName" == "THETA" ];
  then
    echo "aprun -n \${numProcesses} -N \${ppn} --env OMP_NUM_THREADS=\${numOMPthreadsPerRank} -cc depth -d \${numHyperThreadsSkippedPerRank} -j \${numHyperThreadsPerCore} \${@:7:\$#}" >> $SCRATCH/${fileName}/\${scriptName}.sh
  elif [ "$machineName" == "STAMPEDE2" ];
  then
    echo "ibrun \${@:7:\$#}" >> $SCRATCH/${fileName}/\${scriptName}.sh
  elif [ "$machineName" == "PORTER" ];
  then
    if [ "${mpiType}" == "mpi" ];
    then
      mpiexec -n \${numProcesses} \${@:7:\$#}
    elif [ "${mpiType}" == "ampi" ];
    then
      ${BINPATH}charmrun +p1 +vp\${numProcesses} \${@:7:\$#}
    fi
  fi
}


WriteMethodDataForPlotting () {
  for arg in "\${@}"
  do
    echo "echo \"\${arg}\"" >> $SCRATCH/${fileName}/plotInstructions.sh
  done
}

TemporaryDCplotInfo () {
  local scaleRegime=\${1}
  local nodeIndex=0 #\${2} Note that this 6th argument is no longer needed, as these will always print from the original (d,c), regardless of where the local offset is
  local nodeCount=\${3}
  local pDimD=\${4}
  local pDimC=\${5}
  local trickOffset=0 #\${6} Note that this 6th argument is no longer needed, as these will always print from the original (d,c), regardless of where the local offset is
  # New important addition: For special weak scaling, need to print out the number of (d,c) for the binary first, and then each of them in groups of {d,c,(d,c)}
  # Note: still not 100% convinced this is necessary. Need to study scaplot first to make a decision on it.
  # Write to plotInstructions file
  if [ \${scaleRegime} == 2 ];
  then
    echo "echo \"\${nodeCount}\" " >> $SCRATCH/${fileName}/plotInstructions.sh
    curD=\${pDimD}
    curC=\${pDimC}
    trickOffsetTemp=\${trickOffset}
    for ((z=\${nodeIndex}; z<\${nodeCount}; z++))
    do
      echo "echo \"\${curD}\"" >> $SCRATCH/${fileName}/plotInstructions.sh
      echo "echo \"\${curC}\"" >> $SCRATCH/${fileName}/plotInstructions.sh
      echo "echo \"(\${curD},\${curC})\"" >> $SCRATCH/${fileName}/plotInstructions.sh
      trickOffsetTempMod=\$(( trickOffsetTemp % 4 ))
      if [ \${trickOffsetTempMod} == 0 ];
      then
	curD=\$(( \${curD} / 2))
	curC=\$(( \${curC} * 2))
      else
	curD=\$(( \${curD} * 2))
      fi
      trickOffsetTemp=\$(( \${trickOffsetTemp} + 1 ))
    done
  fi
}

WriteMethodDataForCollectingStage1 () {
  local MethodTag=\${1}
  local FileNameBase=\${2}
  local FileName1=\${3}
  local FileName2=\${4}
  local WriteFile=\${5}

  echo "echo \"0\"" >> \${WriteFile}
  echo "echo \"\${MethodTag}\"" >> \${WriteFile}
  echo "echo \"\${FileName1}\"" >> \${WriteFile}
  if [ "\${MethodTag}" != "bscf" ];
  then
    echo "echo \"\${FileName2}\"" >> \${WriteFile}
  fi
  if [ "${profType}" == "PC" ] || [ "${profType}" == "PCT" ];
  then
    echo "echo \"\${FileNameBase}_critter\"" >> \${WriteFile}
  fi
  if [ "${profType}" == "PT" ] || [ "${profType}" == "PCT" ];
  then
    echo "echo \"\${FileNameBase}_timer\"" >> \${WriteFile}
  fi
}

WriteMethodDataForCollectingStage2 () {
  # Because 'Pre' (Stage1) collapses the NumLaunchesPerBinary, we do not want to overcount.
  local launchID=\${1}
  if [ \${launchID} -eq 1 ];
  then
    local MethodTag=\${2}
    local PreFileNameBase=\${3}
    local PreFileName1=\${4}
    local PreFileName2=\${5}
    local PostFileNameBase=\${6}
    local PostFileName1=\${7}
    local PostFileName2=\${8}
    local WriteFile=\${9}

    echo "echo \"0\"" >> \${WriteFile}
    echo "echo \"\${MethodTag}\"" >> \${WriteFile}
    echo "echo \"\${PostFileName1}\"" >> \${WriteFile}
    if [ "\${MethodTag}" != "bscf" ];
    then
      echo "echo \"\${PostFileName2}\"" >> \${WriteFile}
    fi
    if [ "${profType}" == "PC" ] || [ "${profType}" == "PCT" ];
    then
      echo "echo \"\${PostFileNameBase}_critter\"" >> \${WriteFile}
    fi
    if [ "${profType}" == "PT" ] || [ "${profType}" == "PCT" ];
    then
      echo "echo \"\${PostFileNameBase}_timer\"" >> \${WriteFile}
    fi
    echo "echo \"\${PreFileName1}\"" >> \${WriteFile}
    if [ "\${MethodTag}" != "bscf" ];
    then
      echo "echo \"\${PreFileName2}\"" >> \${WriteFile}
    fi
    if [ "${profType}" == "PC" ] || [ "${profType}" == "PCT" ];
    then
      echo "echo \"\${PreFileNameBase}_critter\"" >> \${WriteFile}
    fi
    if [ "${profType}" == "PT" ] || [ "${profType}" == "PCT" ];
    then
      echo "echo \"\${PreFileNameBase}_timer\"" >> \${WriteFile}
    fi
  fi
}

launchJobsPortal () {
  # Launch performance job always.
  launchJobs \${@:2:6} \${1}_PERFORMANCE \${@:8:\${#}}

  # If analysis is turned on, launch Profiling job and Critter job.
  if [ "${profType}" == "PC" ] || [ "${profType}" == "PCT" ];
  then
    launchJobs \${@:2:6} \${1}_CRITTER \${@:8:\${#}}
  fi
  if [ "${profType}" == "PT" ] || [ "${profType}" == "PCT" ];
  then
    launchJobs \${@:2:6} \${1}_TIMER \${@:8:\${#}}
  fi
}

###################################################### Method Definitions ######################################################

# For CA-CQR2
collectPlotTags=()
source Libraries/camfs/cacqr2.sh
source Libraries/camfs/cfr3d.sh
source Libraries/camfs/mm3d.sh
source Libraries/candmc/bsqr.sh
source Libraries/candmc/bscf.sh

###################################################### Method Launches ######################################################

# Note: in future, I may want to decouple numBinaries and numPlotTargets, but only when I find it necessary
# Write to Plot Instructions file, for use by SCAPLOT makefile generator
echo "echo \"1\"" > $SCRATCH/${fileName}/plotInstructions.sh
echo "echo \"${fileNameToProcess}\"" >> $SCRATCH/${fileName}/plotInstructions.sh
echo "echo \"${numTests}\"" >> $SCRATCH/${fileName}/plotInstructions.sh
echo "echo \"${machineName}\"" >> $SCRATCH/${fileName}/plotInstructions.sh
echo "echo \"${profType}\"" >> $SCRATCH/${fileName}/plotInstructions.sh
echo "echo \"${nodeScaleFactor}\"" >> $SCRATCH/${fileName}/plotInstructions.sh

WriteHeaderForCollection $SCRATCH/${fileName}/collectInstructionsStage1.sh
WriteHeaderForCollection $SCRATCH/${fileName}/collectInstructionsStage2.sh

for ((i=1; i<=${numTests}; i++))
do
  echo -e "\nTest #\${i}\n"

  # Nodes
  read -p "Enter starting number of nodes for this test: " startNumNodes
  read -p "Enter ending number of nodes for this test: " endNumNodes

  # Now allowing for WS and SS of any kind in a single bench job.
  read -p "Enter Scaling regime:\
	   [0 -> Weak scaling with increasingly rectangular matrix/grid for QR, regular scaling scheme for CF,MM\
	    1 -> Strong scaling with increasingly rectangular grid for QR, larger cubic grid for CF,MM\
            2 -> Weak scaling with alternating scaling scheme for QR only]: " scaleRegime

  scale="WS"
  if [ \${scaleRegime} == "1" ];		# The rest are WS, which is it already set as
  then
    scale="SS"
  fi

  echo "echo \"\${scale}\"" >> $SCRATCH/${fileName}/plotInstructions.sh

  nodeCount=\$(findCountLength \${startNumNodes} \${endNumNodes} 3 ${nodeScaleFactor})
  echo "echo \"\${nodeCount}\" " >> $SCRATCH/${fileName}/plotInstructions.sh

  curNumNodes=\${startNumNodes}
  for ((j=0; j<\${nodeCount}; j++))
  do
    echo "echo \"\${curNumNodes}\"" >> $SCRATCH/${fileName}/plotInstructions.sh
    curNumNodes=\$(( \${curNumNodes} * ${nodeScaleFactor} ))
  done

  read -p "Enter matrix dimension m: " matrixDimM
  read -p "Enter matrix dimension n: " matrixDimN
  read -p "Enter number of iterations (per launch): " numIterations

  echo "echo \"\${matrixDimM}\"" >> $SCRATCH/${fileName}/plotInstructions.sh
  echo "echo \"\${matrixDimN}\"" >> $SCRATCH/${fileName}/plotInstructions.sh

  j=1
  while [ 1 -eq 1 ];		# Loop iterates until user says stop
  do
    echo -e "\nStage #\${j}"

    # Echo for SCAPLOT makefile generator
    read -p "Enter binary tag [0 for CA-CQR2, 1 for bsqr, 2 for CFR3D, 3 for bscf, 4 for quit, 5 for rsqr, 6 for mm3d]: " binaryTagChoice
    echo "echo \"\${binaryTagChoice}\"" >> $SCRATCH/${fileName}/collectInstructionsStage1.sh
    echo "echo \"\${binaryTagChoice}\"" >> $SCRATCH/${fileName}/collectInstructionsStage2.sh
    echo "echo \"\${binaryTagChoice}\"" >> $SCRATCH/${fileName}/plotInstructions.sh
    # break case
    if [ \${binaryTagChoice} -eq "4" ];
    then
      echo "done with iteration \${i} of ${numTests}"
      break
    fi

    binaryTag=""
    if [ \${binaryTagChoice} == 0 ];
    then
      binaryTag=cqr2
    elif [ \${binaryTagChoice} == 1 ];
    then
      binaryTag=bsqr
    elif [ \${binaryTagChoice} == 2 ];
    then
      binaryTag=cfr3d
    elif [ \${binaryTagChoice} == 3 ];
    then
      binaryTag=bscf
    elif [ \${binaryTagChoice} == 5 ];
    then
      binaryTag=rsqr
    elif [ \${binaryTagChoice} == 6 ];
    then
      binaryTag=mm3d
    fi

    binaryPath=${BINPATH}\${binaryTag}_${machineName}
    if [ "${machineName}" == "PORTER" ];
    then
      binaryPath=\${binaryPath}_${mpiType}
    elif [ "${machineName}" == "BLUEWATERS" ];
    then
      # special case, only for CAMFS, not for bench_scalapack routines
      if [ "\${binaryTag}" == "cqr2" ] || [ "\${binaryTag}" == "cfr3d" ];
      then
        binaryPath=${BINPATH}\${binaryTag}_${machineName}_${GPU}
      fi
    fi

    # State variables that scale with nodes must be initialized here. Otherwise, repeated input is needed in loops below
    if [ \${binaryTag} == 'cqr2' ];
    then
      read -p "Enter start range of starting tunable processor grid dimension c: " startStartPdimC
      read -p "Enter end range of starting tunable processor grid dimension c (for any node count * ppn pairing): " endStartPdimC
      # invCutOff shouldn't be asked for. It should, for now, range up to 2 from 0, unless I am seeing a pattern in performance.
      read -p "Do you want to vary the invCutOff? yes[y] or no[n]: " invCutOffDec
    elif [ \${binaryTag} == 'bsqr' ] || [ \${binaryTag} == 'rsqr' ];
    then
      read -p "Enter the starting number of processor grid columns: " startStartNumPcols
      read -p "Enter the ending number of processor grid columns: " endStartNumPcols
      read -p "Enter the minimum block size: " minBlockSize
      read -p "Enter the maximum block size: " maxBlockSize
      # Anything else? Is above, sufficient?
    elif [ \${binaryTag} == 'cfr3d' ];
    then
      read -p "Enter the starting (cubic) processor grid dimension: " cubeDim
    elif [ \${binaryTag} == 'bscf' ];
    then
      read -p "Enter the minimum block size: " minBlockSize
      read -p "Enter the maximum block size: " maxBlockSize
      # Anything else? Is above, sufficient?
    elif [ \${binaryTag} == 'mm3d' ];
    then
      read -p "Gemm[0] or TRMM[1]: " gemmORtrmmChoice
      read -p "Bcast+Allreduce[0] or Allgather+Allreduce[1]: " bcastORallgatherChoice
      read -p "Enter matrix dimension k: " matrixDimK
      read -p "Enter the starting (cubic) processor grid dimension: " cubeDim
    fi

    for ((curLaunchID=1; curLaunchID<=${NumLaunchesPerBinary}; curLaunchID+=1));
    do
      # Initialize all possible variables that change with node count
      # shared
      nodeIndex=0
      curMatrixDimM=\${matrixDimM}
      curMatrixDimN=\${matrixDimN}
      if [ \${binaryTag} == 'mm3d' ];
      then
        curMatrixDimK=\${matrixDimK}
      fi
      WShelpcounter=0			# change if we want to start at node offset (rare)
      # cqr2
      pDimCArray=()
      pDimCArrayOrig=()
      rangePdimClen=0
      if [ \${binaryTag} == 'cqr2' ];
      then
        for ((w=\${startStartPdimC}; w<=\${endStartPdimC}; w*=2));
        do
          pDimCArray+=(\${w})
          pDimCArrayOrig+=(\${w})
          rangePdimClen=\$(( \${rangePdimClen} + 1 ))
        done
      fi
      # bsqr/rsqr
      numPcolsArray=()
      numPcolsArrayOrig=()
      rangeNumPcolslen=0
      if [ \${binaryTag} == 'bsqr' ] || [ \${binaryTag} == 'rsqr' ];
      then
        for ((w=\${startStartNumPcols}; w<=\${endStartNumPcols}; w*=2));
        do
          numPcolsArray+=(\${w})
          numPcolsArrayOrig+=(\${w})
          rangeNumPcolslen=\$(( \${rangeNumPcolslen} + 1 ))
        done
      fi
      # cfr3d
      curCubeDim=\${cubeDim}
      for ((curNumNodes=\${startNumNodes}; curNumNodes<=\${endNumNodes}; curNumNodes*=${nodeScaleFactor}));
      do
        minPPN=\${ppnMinListRunTime[\${nodeIndex}]}
        maxPPN=\${ppnMaxListRunTime[\${nodeIndex}]}
        for ((curPPN=\${minPPN}; curPPN<=\${maxPPN}; curPPN*=${ppnScaleFactor}));
        do
          numProcesses=\$(( \${curNumNodes} * \${curPPN} ))
          StartingNumProcesses=\$(( \${startNumNodes} * \${curPPN} ))

          minTPR=\${tprMinListRunTime[\${nodeIndex}]}
          maxTPR=\${tprMaxListRunTime[\${nodeIndex}]}
	  for ((curTPR=\${minTPR}; curTPR<=\${maxTPR}; curTPR*=${tprScaleFactor}));
	  do
            # Make sure we are in a suitable range
            numPEsPerNode=\$(( \${curPPN} * \${curTPR} ))
            if [ ${minPEcountPerNode} -le \${numPEsPerNode} ] && [ ${maxPEcountPerNode} -ge \${numPEsPerNode} ];
            then
	      # Now decide on a method:

	      if [ \${binaryTag} == 'cqr2' ];
	      then
		# Below: note that the STARTING dimC is being changed. The parameters that aren't solely dependent on the node count are
		#   changed here and not in launchTag***
		for ((w=0; w<\${rangePdimClen}; w+=1));
		do
		  pDimC=\${pDimCArray[\${w}]}
		  pDimCsquared=\$(( \${pDimC} * \${pDimC} ))
		  pDimD=\$(( \${numProcesses} / \${pDimCsquared} ))

		  # Special check because performance for 16 PPN, 4 TPR shows superior performance for the skinniest grid
                  isSpecial=1

		  # Check if pDimC is too big. If so, pDimD will be 0
		  if [ \${pDimD} -ge \${pDimC} ] && [ \${isSpecial} == 1 ];
		  then
		    originalPdimC=\${pDimCArrayOrig[\${w}]}
		    originalPdimCsquared=\$(( \${originalPdimC} * \${originalPdimC} ))
		    originalPdDimD=\$(( \${StartingNumProcesses} / \${originalPdimCsquared} ))
		    launch\${binaryTag} \${scale} \${binaryPath} \${numIterations} \${curLaunchID} \${curNumNodes} \${curPPN} \${curTPR} \${curMatrixDimM} \${curMatrixDimN} \${matrixDimM} \${matrixDimN} \${originalPdDimD} \${originalPdimC} \${pDimD} \${pDimC} \${nodeIndex} \${scaleRegime} \${nodeCount} \${WShelpcounter} \${invCutOffDec}
		  fi
		done
	      elif [ \${binaryTag} == 'bsqr' ] || [ \${binaryTag} == 'rsqr' ];
	      then
                # Special case to watch out for.
		for ((w=0; w<\${rangeNumPcolslen}; w+=1));
		do
		  numPcols=\${numPcolsArray[\${w}]}
		  numProws=\$(( \${numProcesses} / \${numPcols} ))

                  isSpecial=1

		  if [ \${numPcols} -le \${numProws} ] && [ \${isSpecial} == 1 ];
		  then
		    originalNumPcols=\${numPcolsArrayOrig[\${w}]}
		    originalNumProws=\$(( \${StartingNumProcesses} / \${originalNumPcols} ))
                    sharedBinaryTag="bsqr"	# Even if rsqr, use bsqr and then have the corresponding method use the new argument for binaryTag
		    launch\${sharedBinaryTag} \${binaryTag} \${scale} \${binaryPath} \${numIterations} \${curLaunchID} \${curNumNodes} \${curPPN} \${curTPR} \${curMatrixDimM} \${curMatrixDimN} \${matrixDimM} \${matrixDimN} \${originalNumProws} \${originalNumPcols} \${numProws} \${minBlockSize} \${maxBlockSize} \${nodeIndex} \${scaleRegime} \${nodeCount}
		  fi
		done
	      elif [ \${binaryTag} == 'cfr3d' ];
	      then
		launch\${binaryTag} \${scale} \${binaryPath} \${numIterations} \${curLaunchID} \${curNumNodes} \${curPPN} \${curTPR} \${curMatrixDimM} \${matrixDimM} \${cubeDim} \${curCubeDim} \${nodeIndex} \${scaleRegime} \${nodeCount}
	      elif [ \${binaryTag} == 'bscf' ];
	      then
		launch\${binaryTag} \${scale} \${binaryPath} \${numIterations} \${curLaunchID} \${curNumNodes} \${curPPN} \${curTPR} \${curMatrixDimM} \${matrixDimM} \${minBlockSize} \${maxBlockSize} \${nodeIndex} \${scaleRegime} \${nodeCount}
	      elif [ \${binaryTag} == 'mm3d' ];
	      then
		launch\${binaryTag} \${scale} \${binaryPath} \${numIterations} \${curLaunchID} \${curNumNodes} \${curPPN} \${curTPR} \${gemmORtrmmChoice} \${bcastORallgatherChoice} \${curMatrixDimM} \${curMatrixDimN} \${curMatrixDimK} \${matrixDimM} \${matrixDimN} \${matrixDimK} \${cubeDim} \${curCubeDim} \${nodeIndex} \${scaleRegime} \${nodeCount}
	      fi
            fi
          done
        done
        # Update all loop variables for each increasing node count
        # Assuming that pDimD can always just be calculated from pDimC and NumProcesses as we scale
        if [ \${scaleRegime} == 0 ];
        then
	  # below: shared
	  curMatrixDimM=\$(( \${curMatrixDimM} * 2 ))
	  # below: ca-cqr2
	  if [ \${binaryTag} == 'cqr2' ];
          then
	    #pDimD=\$(( \${pDimD} * 2 ))
	    echo "Do nothing"
          fi
          # below: bench scala qr
	  if [ \${binaryTag} == 'bsqr' ] || [ \${binaryTag} == 'rsqr' ];
          then
	    echo "Do nothing"
	  fi
	  # below: cfr3d
	  if [ \${binaryTag} == 'cfr3d' ] || [ \${binaryTag} == 'mm3d' ];
          then
            curCubeDim=\$(( \${curCubeDim} * 2 ))
            if [ \${binaryTag} == 'mm3d' ];
	    then
	      curMatrixDimN=\$(( \${curMatrixDimN} * 2 ))
	      curMatrixDimK=\$(( \${curMatrixDimK} * 2 ))
	    fi
          fi
        elif [ \${scaleRegime} == 1 ];
	then
	  # below: ca-cqr2
	  if [ \${binaryTag} == 'cqr2' ];
          then
            #pDimD=\$(( \${pDimD} * 2 ))
	    echo "Do nothing"
	  fi
	  # below: bench scala qr
	  if [ \${binaryTag} == 'bsqr' ] || [ \${binaryTag} == 'rsqr' ];
          then
	    echo "Do nothing"
	  fi
	  # below: cfr3d
	  if [ \${binaryTag} == 'cfr3d' ] || [ \${binaryTag} == 'mm3d' ];
          then
            curCubeDim=\$(( \${curCubeDim} * 2 ))
          fi
        elif [ \${scaleRegime} == 2 ];
	then
	  immWS=\$(( \${WShelpcounter} % 4 ))
	  if [ \${immWS} == 0 ];
	  then
            # shared
	    curMatrixDimM=\$(( \${curMatrixDimM} / 2 ))
	    curMatrixDimN=\$(( \${curMatrixDimN} * 2 ))
	    # below: ca-cqr2
	    if [ \${binaryTag} == 'cqr2' ];
            then
	      #pDimD=\$(( \${pDimD} / 2 ))
	      for ((w=0; w<\${rangePdimClen}; w+=1));
	      do
	        pDimC=\${pDimCArray[\${w}]}
	        pDimCArray[\${w}]=\$(( \${pDimC} * 2 ))	# update
	      done
            fi
            # below: bench scala qr
	    if [ \${binaryTag} == 'bsqr' ] || [ \${binaryTag} == 'rsqr' ];
            then
              for ((w=0; w<\${rangeNumPcolslen}; w+=1));
              do
                numPcols=\${numPcolsArray[\${w}]}
                numPcolsArray[\${w}]=\$(( \${numPcols} * 2 ))
              done
            fi
	  else
            # shared
	    curMatrixDimM=\$(( \${curMatrixDimM} * 2 ))
	    # below: ca-cqr2
	    if [ \${binaryTag} == 'cqr2' ];
            then
              #pDimD=\$(( \${pDimD} * 2 ))
	      echo "Do nothing"
	    fi
            # below: bench scala qr
	    if [ \${binaryTag} == 'bsqr' ] || [ \${binaryTag} == 'rsqr' ];
            then
	      echo "Do nothing"
	    fi
          fi
          WShelpcounter=\$(( \${WShelpcounter} + 1 ))
	fi
	nodeIndex=\$(( \${nodeIndex} + 1 ))
      done
    done
    j=\$(( \${j} + 1 ))
    echo "echo \"1\"" >> $SCRATCH/${fileName}/collectInstructionsStage1.sh	# Signals end of the data files for this specific methodID
    echo "echo \"1\"" >> $SCRATCH/${fileName}/collectInstructionsStage2.sh	# Signals end of the data files for this specific methodID
    echo "echo \"1\"" >> $SCRATCH/${fileName}/plotInstructions.sh	# Signals end of the data files for this specific methodID
  done
done
EOF


bash $SCRATCH/${fileName}.sh
#rm $SCRATCH/${fileName}.sh

# Copy a local version to Scripts directory so that it can be used on the local side to generate plots.
# But its important that we keep a backup in SCRATCH/fileName in case we overwrite collectInstructionsStage1.sh, we can always write it back.
# cp $SCRATCH/${fileName}/collectInstructionsStage1.sh collectInstructionsStage1.sh		// collectInstructionsStage1.sh will not be needed locally anymore.
# Do not copy collectInstructionsStage2 to local directory.

# Note that for Porter, no need to do this, since we are submitting to a queue
if [ "${machineName}" == "BGQ" ] || [ "${machineName}" == "THETA" ] || [ "${machineName}" == "STAMPEDE2" ] || [ "${machineName}" == "BLUEWATERS" ];
then
  mkdir $SCRATCH/${fileName}/bin
  mv ../bin/* $SCRATCH/${fileName}/bin
  cd $SCRATCH

  # Submit all scripts
  curLaunchID=1
  while [ ${curLaunchID} -le ${NumLaunchesPerBinary} ];
  do
    curNumNodes=${minNumNodes}
    ppnIndex=0
    while [ ${curNumNodes} -le ${maxNumNodes} ];
    do
      minPPN=${ppnMinList[${ppnIndex}]}
      maxPPN=${ppnMaxList[${ppnIndex}]}
      curPPN=${minPPN}
      tprIndex=0
      while [ ${curPPN} -le ${maxPPN} ];
      do
        minTPR=${tprMinList[${tprIndex}]}
        maxTPR=${tprMaxList[${tprIndex}]}
        curTPR=${minTPR}
        while [ ${curTPR} -le ${maxTPR} ];
        do
          # Make sure we are in a suitable range
          numPEsPerNode=$(( ${curPPN} * ${curTPR} ))
          if [ ${minPEcountPerNode} -le ${numPEsPerNode} ] && [ ${maxPEcountPerNode} -ge ${numPEsPerNode} ];
          then
            if [ "${machineName}" == "BGQ" ] || [ "${machineName}" == "THETA" ];
            then
              qsub ${fileName}/script_${fileID}id_${roundID}round_${curLaunchID}launchID_${curNumNodes}nodes_${curPPN}ppn_${curTPR}tpr.sh
            elif [ "${machineName}" == "BLUEWATERS" ];
            then
              qsub ${fileName}/script_${fileID}id_${roundID}round_${curLaunchID}launchID_${curNumNodes}nodes_${curPPN}ppn_${curTPR}tpr.pbs
            else
              chmod +x ${fileName}/script_${fileID}id_${roundID}round_${curLaunchID}launchID_${curNumNodes}nodes_${curPPN}ppn_${curTPR}tpr.sh
              sbatch --mail-user=${MyEmail} --mail-type=all ${fileName}/script_${fileID}id_${roundID}round_${curLaunchID}launchID_${curNumNodes}nodes_${curPPN}ppn_${curTPR}tpr.sh
            fi
          fi
          curTPR=$(( ${curTPR} * ${ppnScaleFactor} ))
        done
        curPPN=$(( ${curPPN} * ${tprScaleFactor} ))
        tprIndex=$(( ${tprIndex} + 1 ))
      done
      curNumNodes=$(( ${curNumNodes} * ${nodeScaleFactor} ))
      ppnIndex=$(( ${ppnIndex} + 1 ))
    done
    curLaunchID=$(( ${curLaunchID} + 1 ))
  done
fi
