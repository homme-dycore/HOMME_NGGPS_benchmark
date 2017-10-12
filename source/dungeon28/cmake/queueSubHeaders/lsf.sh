#!/bin/bash

#BSUB -a poe
#BSUB -P ${HOMME_PROJID}

#BSUB -q premium
#BSUB -W 0:15
#BSUB -x

#BSUB -J ${TEST_NAME}

#BSUB -o ${TEST_NAME}.stdout.%J
#BSUB -e ${TEST_NAME}.stderr.%J

#BSUB -n ${NUM_CPUS}
#BSUB -R "span[ptile=${NUM_CPUS}]" 

