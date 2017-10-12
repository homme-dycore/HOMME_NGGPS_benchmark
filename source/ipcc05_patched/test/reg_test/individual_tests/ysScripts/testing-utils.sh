setTestDirs() {

  # Determine some locations
  DIR_NAME=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd -P)

  # Determine the location of the tests results (for yellowstone)
  RESULT_DIR=$(cd "${DIR_NAME}/../../results/yellowstone" && pwd -P)

  # Set the location of the "build" base directory
  if [ -n "$1" -a -d "$1" ]; then
    # Use this location as the base of the file structure for the tests
    BUILD_DIR=$1
  else
    # Set the build directory from the set file structure
    BUILD_DIR=$(cd `dirname $DIR_NAME/../../../../..` && pwd -P)/build
  fi

}

yellowstoneLSFFile() {

  RUN_SCRIPT=$1

  #delete the file if it exists
  rm -f $RUN_SCRIPT

  # Set up some yellowstone boiler plate
  echo "#!/bin/bash" >> $RUN_SCRIPT
  echo ""  >> $RUN_SCRIPT # newlines

  echo "#BSUB -a poe" >> $RUN_SCRIPT

  # To do: move this check up and properly handle the error status
  if [ -n "$HOMME_PROJID" ]; then
    echo "#BSUB -P $HOMME_PROJID" >> $RUN_SCRIPT
  else
    echo "PROJECT CHARGE ID (HOMME_PROJID) not set"
    exit -1
  fi 

  echo "" >> $RUN_SCRIPT # newline

  echo "#BSUB -q small" >> $RUN_SCRIPT
  echo "#BSUB -W 0:20" >> $RUN_SCRIPT

  echo "" >> $RUN_SCRIPT # newline

  echo "#BSUB -x" >> $RUN_SCRIPT

  echo "" >> $RUN_SCRIPT # newline

  echo "#BSUB -R \"select[scratch_ok > 0 ]\"" >> $RUN_SCRIPT

  echo "" >> $RUN_SCRIPT # newline

  # Set the job name
  echo "#BSUB -J $testName" >> $RUN_SCRIPT
  echo "" >> $RUN_SCRIPT

  # Set the output and error filenames
  echo "#BSUB -o $testName.stdout.%J" >> $RUN_SCRIPT
  echo "#BSUB -e $testName.stderr.%J" >> $RUN_SCRIPT
  echo "" >> $RUN_SCRIPT

  # Set the ncpus and ranks per MPI
  echo "#BSUB -n $num_cpus" >> $RUN_SCRIPT
  echo '#BSUB -R "span[ptile='$num_cpus']" ' >> $RUN_SCRIPT

  echo "" >> $RUN_SCRIPT

  echo "cd $outputDir" >> $RUN_SCRIPT

  echo "" >> $RUN_SCRIPT

}


yellowstoneExec() {
  RUN_SCRIPT=$1
  EXEC=$2
  echo "mpirun.lsf $EXEC" >> $RUN_SCRIPT

}

printSubmissionSummary() {

  # Output a summary of the test name along with the bsub job id for easy reference
  echo "" # newline
  echo "############################################################################"
  echo "Summary of submissions"
  echo "############################################################################"
  for i in $(seq 0 $(( ${#SUBMIT_TEST[@]} - 1)))
  do
    echo "Test ${SUBMIT_TEST[$i]} has ID ${SUBMIT_JOB_ID[i]}"
  done
  echo "############################################################################"

}

queueWait() {

  for i in $(seq 0 $(( ${#SUBMIT_TEST[@]} - 1)))
  do
    echo -n "Examining status of ${SUBMIT_TEST[$i]}..."
    jobID=${SUBMIT_JOB_ID[i]}
    jobFinished=false

    while ! $jobFinished;
    do
      # Test if the job exists
      jobStat=`bjobs -a $jobID | tail -n 1 | awk '{print $3}'`

      # Print the status of the job
      echo -n "$jobStat..."

      # if the job is registered in the queue and the status is PEND or RUN then wait
      if [ -n "$jobStat" -a "$jobStat" == "PEND" -o "$jobStat" == "RUN" ]; then
        # Job still in queue or running
        sleep 60 # sleep for 60s
      else # if jobStat=DONE, EXIT or it is finished and no longer registered in the queue
        jobFinished=true
        echo "FINISHED..."
      fi
    done
  done
  echo "############################################################################"

}

diffStdOut() {

  # Should be a unique file
  diffFile="diff.${SUBMIT_JOB_ID[0]}"
  echo "Concatenating all diff output into $diffFile"

  # Then diff with the stored results (yellowstone only)
  for i in $(seq 0 $(( ${#SUBMIT_TEST[@]} - 1)))
  do
    THIS_TEST=${SUBMIT_TEST[$i]}
    # The following is not very clean
    NEW_RESULT=${DIR_NAME}/${THIS_TEST}.stdout.${SUBMIT_JOB_ID[$i]}
    SAVED_RESULT=${RESULT_DIR}/${THIS_TEST}/${THIS_TEST}.stdout
    #stripAppendage $NEW_RESULT
    echo "diff $NEW_RESULT $SAVED_RESULT" >> $diffFile
    # append the output to 
    diff $NEW_RESULT $SAVED_RESULT >> $diffFile
  done
}

submitTestsToLSF() {

  echo "Submitting ${num_submissions} jobs to queue"

  SUBMIT_TEST=()
  SUBMIT_JOB_ID=()
  SUBMIT_TEST=()

  # Loop through all of the tests
  for subNum in $(seq 1 ${num_submissions})
  do

    subFile=subFile${subNum}
    subFile=${!subFile}
    echo "subFile=${subFile}"
    subJobName=`basename ${subFile} .sh`
    echo "subJobName=$subJobName"

    # setup file for stdout and stderr redirection
    THIS_STDOUT=${subJobName}.out
    THIS_STDERR=${subJobName}.err

    # Run the command
    # For some reason bsub must not be part of a string
    echo -n "Submitting test ${subJobName} to the queue... "
    bsub < ${subFile} > $THIS_STDOUT 2> $THIS_STDERR
    BSUB_STAT=$?

    # Do some error checking
    if [ $BSUB_STAT == 0 ]; then
      # the command was succesful
      BSUB_ID=`cat $THIS_STDOUT | awk '{print $2}' | sed  's/<//' | sed -e 's/>//'`
      echo "successful job id = $BSUB_ID"
      SUBMIT_TEST+=( "${subJobName}" )
      SUBMIT_JOB_ID+=( "$BSUB_ID" )
    else 
      echo "failed with message:"
      cat $THIS_STDERR
      exit -1
    fi
    rm $THIS_STDOUT
    rm $THIS_STDERR
  done
}

parseStdout() {

  # Then diff with the stored results (yellowstone only)
  for i in $(seq 0 $(( ${#SUBMIT_TEST[@]} - 1)))
  do
    THIS_TEST=${SUBMIT_TEST[$i]}

    # Need to remove "-run" from the test name
    #   This is an ugly hack but otherwise this takes a lot of reformatting
    THIS_TEST=`echo $THIS_TEST | sed 's/-run//'`

    # The following is not very clean
    NEW_RESULT=${HOMME_TESTING_DIR}/${THIS_TEST}.stdout.${SUBMIT_JOB_ID[$i]}
    PARSE_RESULT=${HOMME_TESTING_DIR}/${THIS_TEST}.stdout

    if [ "$HOMME_Submission_Type" = lsf ]; then
      stripAppendage $NEW_RESULT
    fi
    grep -e '=' ${NEW_RESULT} | grep -iv bsub > ${PARSE_RESULT}

  done
}
