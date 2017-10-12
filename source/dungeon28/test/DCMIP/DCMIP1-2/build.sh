#!/bin/tcsh
set TEST_DIR = `pwd`
cd ../../..
make -j6 preqx_L64
cd $TEST_DIR
