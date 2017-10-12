#!/bin/bash

# The lists of tests to run
source lsf-list.sh

# The testing utilities
source testing-utils.sh

# Submit the tests to the queue
submitTestsToLSF

# Print a summary of the submissions
printSubmissionSummary

# Wait for the jobs to run through the queue
queueWait

# Diff the output files with those saved in the repo
diffStdOut


