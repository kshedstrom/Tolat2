#  ******************************************************************
#  Kate Hedstrom  (kate@ahab.rutgers.edu)
#  ******************************************************************
#
# I always use -g here because I often end up running grid
# under a debugger.  This may not be true for you, but be
# VERY careful about optimizing the grid program.
#
# CPPFLAGS      Extra flags for the C preprocessor
# CPP           Name of the C preprocessor
# LDR           Program to load the objects into an executable
# LDFLAGS       Flags to the loader
# RANLIB        Name of ranlib command (SysV doesn't have)
# LIBDIR        Directory where libspem.a goes
# CDFLIB        Name of netCDF library, either -lnetcdf or null
#
   NETCDF_ROOT = /usr/local/pkg/MPI/GCC/11.3.0/OpenMPI/4.1.4/netCDF-C++4/4.3.1
  NETCDF_ROOT2 = /usr/local/pkg/MPI/GCC/11.3.0/OpenMPI/4.1.4/netCDF/4.9.0
      CPPFLAGS = -I$(NETCDF_ROOT)/include
        LIBDIR = $(NETCDF_ROOT)/lib
       LIBDIR2 = $(NETCDF_ROOT2)/lib
      CXXFLAGS = -g
        RANLIB = ranlib

       LDFLAGS =
        CDFLIB = -L$(LIBDIR) -lnetcdf_c++4 -L$(LIBDIR2) -lnetcdf

Tolat2: Proj.o Tolat2.o
	$(CXX) -o Tolat2 $(LDFLAGS) $(CXXFLAGS) Proj.o Tolat2.o $(CDFLIB)

Coast: Proj.o Coast.o
	$(CXX) -o Coast $(LDFLAGS) $(CXXFLAGS) Proj.o Coast.o

Coast_2: Proj.o Coast_2.o
	$(CXX) -o Coast_2 $(LDFLAGS) $(CXXFLAGS) Proj.o Coast_2.o

Toxy: Toxy.o
	$(CXX) -o Toxy $(LDFLAGS) $(CXXFLAGS) Toxy.o -lnetcdf_g++ -lnetcdf_gcc

Toxy2: Toxy2.o
	$(CXX) -o Toxy2 $(LDFLAGS) $(CXXFLAGS) Toxy2.o -lnetcdf_g++ -lnetcdf_gcc

depend:
	g++ -M *.cc > makedepend

clean:
	rm -f *.o gmeta core *.trace

include makedepend
