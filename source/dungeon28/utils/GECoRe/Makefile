# NETCDF library directories
NETCDF_INCLUDEDIR=/usr/include
NETCDF_LIBDIR=/usr/lib

# Turn on BLAS support
USEBLAS = True

# Library files to include
LDFILES= -lnetcdf -lnetcdf_c++

# Local files
FILES= gecore.cpp \
		Interpolator.cpp Interpolator_search.cpp HOMMEInterpolator.cpp \
		LatitudeLongitudeGrid.cpp CubedSphereGrid.cpp HOMMEGrid.cpp \
		Announce.cpp CubedSphereTrans.cpp \
        Preferences.cpp LinearAlgebra.cpp GaussLobattoQuadrature.cpp

# Load system-specific defaults
include Make.defs

CFLAGS+= -I$(NETCDF_INCLUDEDIR)
LDFLAGS+= -L$(NETCDF_LIBDIR)

##
## Build instructions
##
gecore: $(FILES:%.cpp=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

##
## Clean
##
clean:
	rm -f gecore *.o
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

##
## Include dependencies
##
include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE

