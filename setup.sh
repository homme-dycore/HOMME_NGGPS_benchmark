#!/bin/bash

# This script will look for a file called '<system>.config' in the systems directory, and load
# settings if found and quitting if not.  To add a new system, just create a new file.  Currently
# the Cori (NERSC), Cheyenne (NCAR) and Yellowstone (NCAR) systems are supported
SYSTEM="yellowstone"

# The cases we want to run, and whether to submit automatically on running this script:
CASES="small, medium, large"
RUN="true"



main() {
  CheckArguments $1
  LoadConfiguration

  BuildExe original  
  BuildExe optimized 

  MakeCases
  RunCases
}

###################################################################
# CheckArguments Function - verifies we were called correctly and #
#                           that the directory doesn't yet exist, #
#                           then sets source / target dir vars.   #
###################################################################
CheckArguments() {
if [ "$#" -ne 1 ]; then
  echo "Use: ./setup <directory> "
  exit 1
fi

if [ -d $1 ]; then
  echo "Error: Directory '$1' already exists; exiting."
  exit 1
fi

# Set up our source and target directory variables:
SOURCEDIR=$( dirname "$(readlink -f "$0")" )    # Directory of this script
TARGETDIR=$( readlink -f $1 )                   # Benchmark directory

# Remove commas from cases
CASES=${CASES//,/ }
}


###################################################################
# LoadConfiguration Function - Reads a configuration from the     #
#                              systems directory, failing if it   #
#                              doesn't exist.                     #
###################################################################
LoadConfiguration() {
if [ -f "${SOURCEDIR}/systems/${SYSTEM}.config" ]; then
  echo "Using systems/${SYSTEM}.config..."
  source ${SOURCEDIR}/systems/${SYSTEM}.config
else
  echo "Error: No configuration file found for system '${SYSTEM}'"
  exit 1
fi
}


###################################################################
# BuildExes Function - Makes directories for builds, then sets    #
#                      up and executes a build script for both    #
#                      versions of the code.                      #
###################################################################
BuildExe() {
VERSION=$1

mkdir -p ${TARGETDIR}/builds/${VERSION}
local HOMME_ROOT=${SOURCEDIR}/source/${VERSION}
local CONFIG_FILE=${TARGETDIR}/builds/${VERSION}/config.sh

cat << EOF > ${CONFIG_FILE}
#!/bin/bash

rm -rf CMakeFiles CMakeCache.txt

${PRE_RUN}

export HOMME_ROOT=${HOMME_ROOT}

cmake -DQSIZE_D=10 -DPREQX_PLEV=128 -DPREQX_NP=4 \\
  -DCMAKE_SYSTEM_NAME=${OS} \\
  -DCMAKE_Fortran_COMPILER=${FC} \\
  -DCMAKE_C_COMPILER=${CC} \\
  -DCMAKE_CXX_COMPILER=${CXX} \\
  -DNETCDF_DIR=${NETCDF} \\
  -DHOMME_PROJID=${PROJID} \\
  -DENABLE_NANOTIMERS=TRUE \\
  -DUSE_BIT64=TRUE \\
  -DBUILD_HOMME_SWDGX=FALSE \\
  -DBUILD_HOMME_SWEQX=FALSE \\
  -DBUILD_HOMME_PRIMDGX=FALSE \\
  -DPREQX_USE_ENERGY=FALSE preqx \\
  -DFORCE_Fortran_FLAGS=${FFLAGS} \\
  -DFORCE_C_FLAGS=${CFLAGS} \\
  -DFORCE_LINKER_FLAGS=${LFLAGS} \\
  ${HOMME_ROOT}

make -j 6 preqx

EOF

cd ${TARGETDIR}/builds/${VERSION}
chmod +x ./config.sh
echo "Building '${VERSION}' code..."
./config.sh &> build.log 
if [ ! -f ${TARGETDIR}/builds/${VERSION}/src/preqx/preqx ]; then
  echo "Error building preqx executable; quitting."
  exit 1
fi

}

###################################################################
# MakeRunScript Function - Makes the run script for this system   #
###################################################################
MakeRunScript() {
RUN_FILE=$1

cat << EOF > ${RUN_FILE}
${JOBHEADER}
${PRE_RUN}

${JOBSETTINGS}

cd ${TARGETDIR}/cases/${CASE}

if [ -f /tmp/original.exe ]; then
  EXEDIR="/tmp"
else
  EXEDIR="${TARGETDIR}/cases/${CASE}"
fi

for r in \$(seq -f %02g 1 3) ; do
  ${MPIRUN} \${EXEDIR}/original.exe < ${TARGETDIR}/cases/${CASE}/original.nl  > out.original.\${r}.txt
  mv HommeTime_stats HommeTime_stats.original.\${r}
  sleep 1

  ${MPIRUN} \${EXEDIR}/optimized.exe < ${TARGETDIR}/cases/${CASE}/optimized.nl > out.optimized.\${r}.txt
  mv HommeTime_stats HommeTime_stats.optimized.\${r}
  sleep 1
done

# Calculate relative differences:
${SOURCEDIR}/scripts/rel_diff out.original.01.txt ${SOURCEDIR}/verification/${CASE}_original.txt   > reldiff_original.txt
${SOURCEDIR}/scripts/rel_diff out.optimized.01.txt ${SOURCEDIR}/verification/${CASE}_optimized.txt > reldiff_optimized.txt

# Verify results to within set tolerance:
${SOURCEDIR}/scripts/verify reldiff_original.txt   > results.txt
${SOURCEDIR}/scripts/verify reldiff_optimized.txt  >> results.txt


# Write out timing stats:
${SOURCEDIR}/scripts/timing  >> results.txt

# Save results into main directory:
cp results.txt ${TARGETDIR}/results_${CASE}.txt

EOF

# Switch out the number of nodes and ranks in the script we just made:
sed -i "s|NUMBER_OF_NODES|${CONFIG[1]}|g;s|NUMBER_OF_RANKS|${CONFIG[2]}|g" ${RUN_FILE}

}

###################################################################
# GetCaseConfig Function - Gets the number of nodes and ranks     #
#                          for this case based on its name.       #
###################################################################
GetCaseConfig() {

CONFIG=()
if [ -f ${SOURCEDIR}/inputs/${CASE}.nl ] ; then
   for cfg in ${CONFIGS[@]}; do
     cfg_array=($(echo $cfg | sed 's/,/\n/g'))
     if [ "${cfg_array[0]}" == "${CASE}" ] ; then
        CONFIG=(${cfg_array[@]})
     fi
  done
else
  echo "No input namelists found ( ${SOURCEDIR}/inputs/${CASE}.nl ) ; quitting."
  exit 1
fi

if [ ${CONFIG[0]} == "" ]; then
  echo "No configuration found for '${CASE}' case - please add one in the "
  echo "'${SOURCEDIR}/systems/${SYSTEM}.config' file.  Quitting."
  exit 1
fi
}



###################################################################
# MakeCases Function - Makes the case directories, with the       #
#                      necessary input files and job scripts      #
###################################################################
MakeCases() {
for CASE in ${CASES}; do
  if [ -f ${SOURCEDIR}/inputs/${CASE}.nl ] ; then
    # Get the configuration of this case:
    GetCaseConfig 

    # Make the case and case/movies directories:
    mkdir -p ${TARGETDIR}/cases/${CASE}
    mkdir -p ${TARGETDIR}/cases/${CASE}/movies

    # Copy over the vcoord input files:
    cp -rp ${SOURCEDIR}/inputs/vcoord ${TARGETDIR}/cases/${CASE}
    
    # Create the namelists:
    cp ${SOURCEDIR}/inputs/${CASE}.nl ${TARGETDIR}/cases/${CASE}/optimized.nl
    sed 's|horz_num_threads = 1|NThreads = 1|g;s|vert_num_threads = 1||g;s|tracer_num_threads = 1||g' \
        ${SOURCEDIR}/inputs/${CASE}.nl >  ${TARGETDIR}/cases/${CASE}/original.nl

    # Link the executables:
    ln -s ${TARGETDIR}/builds/original/src/preqx/preqx ${TARGETDIR}/cases/${CASE}/original.exe
    ln -s ${TARGETDIR}/builds/optimized/src/preqx/preqx ${TARGETDIR}/cases/${CASE}/optimized.exe

    # Create the run script:
    MakeRunScript ${TARGETDIR}/cases/${CASE}/run.script

  else
    echo "Input files for '${CASE}' case not found; skipping."
  fi
done
}


###################################################################
# RunCases Function - Run s the case directories, with the       #
#                      necessary input files and job scripts      #
###################################################################
RunCases() {
if [ "${RUN}" == "true" ]; then
  for CASE in ${CASES}; do
    if [ -f ${TARGETDIR}/cases/${CASE}/run.script ]; then
       echo "Submitting : ${TARGETDIR}/cases/${CASE}/run.script"
       eval "${JOBSUBMIT} ${TARGETDIR}/cases/${CASE}/run.script"
    fi
  done
fi
}


# Call main script
main $1


