
# Automatically detect scripts and namelist files

file(GLOB NCL_SCRIPTS    "*.ncl")   # get ncl-scripts
file(GLOB SHELL_SCRIPTS  "*.sh")    # get shell-scripts
file(GLOB NAMELISTS      "*.nl")    # get namelists
file(GLOB PYTHON_SCRIPTS "*.py")    # get python scripts
file(GLOB JOB_SCRIPTS    "*.job")   # get jobscripts

# Use make install to copy files to the build directory

install(PROGRAMS ${NCL_SCRIPTS} ${SHELL_SCRIPTS} ${NAMELISTS} ${PYTHON_SCRIPTS} ${JOB_SCRIPTS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
