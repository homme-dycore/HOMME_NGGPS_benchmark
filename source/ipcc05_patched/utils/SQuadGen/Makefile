##############################################################################
# Compiler and flags
CC= g++
CFLAGS= -O3

# NETCDF library directories
NETCDF_INCLUDEDIR=/usr/local/include
NETCDF_LIBDIR=/usr/local/lib

##############################################################################
# DO NOT MODIFY BELOW THIS LINE
##############################################################################

# Library files to include
LDFILES= -lnetcdf -lnetcdf_c++

# Local files
FILES= HOMMEGridGen.cpp GridElements.cpp CSRefinementMap.cpp \
       CubedSphereGrid.cpp IcosahedralFlagGrid.cpp RefineGrid.cpp \
       RefinementTemplateCUBIT.cpp \
       RefinementTemplateLOWCONN.cpp \
       RefinementTemplateLOWCONNOLD.cpp \
       SpringDynamics.cpp \
       Tessellate.cpp \
       lodepng.cpp

# Load system-specific defaults
CFLAGS+= -I$(NETCDF_INCLUDEDIR)
LDFLAGS+= -L$(NETCDF_LIBDIR)

include Make.defs

##
## Build instructions
##
SQuadGen: $(FILES:%.cpp=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

##
## Clean
##
clean:
	rm -f SQuadGen *.o
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

##
## Include dependencies
##
include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE

