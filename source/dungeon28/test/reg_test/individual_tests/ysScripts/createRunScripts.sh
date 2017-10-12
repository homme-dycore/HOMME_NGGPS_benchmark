#!/bin/bash

# Set for now
HOMME_DIR=/glade/scratch/jamroz/ys/homme/homme-cmake/myBuild

source ${HOMME_DIR}/tests/test_list.sh

source testing-utils.sh

lsfListFile=lsf-list.sh

touch $lsfListFile
echo "num_submissions=$num_test_files" > $lsfListFile

for testFileNum in $(seq 1 $num_test_files)
do

  testFile=test_file${testFileNum}
  source ${!testFile}

  testName=`basename ${!testFile} .sh`

  echo "Test $testName has $num_tests pure MPI tests"
  if [ -n "$omp_num_tests" ]; then
    echo "  and $omp_num_tests Hybrid MPI + OpenMP tests"
  fi

  # Create the run script
  thisRunScript=`dirname ${!testFile}`/$testName-run.sh

  outputDir=`dirname ${!testFile}`

  # Set up some yellowstone stuff
  yellowstoneLSFFile $thisRunScript

  for testNum in $(seq 1 $num_tests)
  do
    testExec=test${testNum}
    echo "# Pure MPI test ${testNum}" >> $thisRunScript
    #echo "mpiexec -n $num_cpus ${!testExec} > $testName.out 2> $testName.err" >> $thisRunScript
    yellowstoneExec $thisRunScript "${!testExec}"
    echo "" >> $thisRunScript # new line
  done

  if [ -n "$omp_num_tests" ]; then
    echo "export OMP_NUM_THREADS=$omp_number_threads" >> $thisRunScript
    echo "" >> $thisRunScript # new line
    for testNum in $(seq 1 $omp_num_tests)
    do
       testExec=omp_test${testNum}
       echo "# Hybrid test ${testNum}" >> $thisRunScript
       #echo "mpiexec -n $omp_num_mpi ${!testExec} > $testName.out 2> $testName.err" >> $thisRunScript
       yellowstoneExec $thisRunScript "${!testExec}"
       echo "" >> $thisRunScript # new line
    done
  fi

  echo "subFile$testFileNum=$thisRunScript" >>  $lsfListFile

  # Reset the variables (in case they are not redefined in the next iteration)
  unset omp_num_tests
  unset num_tests

done
