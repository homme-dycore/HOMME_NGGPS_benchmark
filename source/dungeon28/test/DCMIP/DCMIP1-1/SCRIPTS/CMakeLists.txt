# Automatically detect scripts and namelist files

file(GLOB NCL       "*.ncl")   # get ncl-scripts
file(GLOB SHELL     "*.sh")    # get shell-scripts
file(GLOB NAMELISTS "*.nl")    # get namelists
file(GLOB PYTHON    "*.py")    # get python scripts
file(GLOB JOBS      "*.job")   # get jobscripts
file(GLOB ASCII     "*.ascii") # get ascii files

# Use make install to copy files to the build directory

install(PROGRAMS ${NCL} ${SHELL} ${NAMELISTS} ${PYTHON} ${JOBS} ${ASCII} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
