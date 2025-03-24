if DEBUG
AM_FCFLAGS = -g -fcheck=all -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=all -Wuninitialized -finit-real=snan -DDEBUG_TIMESTEP -DDEBUG_NEGATIVE_DEPTH -DEBUG_NEGATIVE_CONC -DDEBUG_EXCESS_CONC -DDEBUG_SPD
AM_CXXFLAGS = -g
else
AM_FCFLAGS = -O2 -falign-loops=16 -march=native -mtune=native
if LINUX
# These optimisation flags are only supported by gcc on linux.
AM_FCFLAGS += -march=native -mtune=native
endif
AM_CXXFLAGS = -O2 -falign-loops=16 
endif
if PARAL
AM_FCFLAGS += -ftree-vectorize -ftree-parallelize-loops=2 -fopt-info-loop -fopenmp
endif

bin_PROGRAMS = kestrel
include_HEADERS = RasterData.h cUTM.h
kestrel_SOURCES = cUTM.cpp RasterData.cpp SetPrecision.f90 Blur.f90 \
	Messages.f90 varStringClass.f90 utilities.f90 Interp2d.f90 \
	utm.f90 RunSettings.f90 GeoTiffRead.f90 DomainSettings.f90 \
	InitConds.f90 Closures.f90 TopogFuncs.f90 Parameters.f90 Limiters.f90 \
	SolverSettings.f90 OutputSettings.f90 TopogSettings.f90 Input.f90 Grid.f90 \
    Equations.f90 HydraulicRHS.f90 MorphodynamicRHS.f90 Redistribute.f90 \
    dem.f90 UpdateTiles.f90 SetSources.f90 NetCDFUtils.f90 \
    Output.f90 Restart.f90 TimeStepper.f90 cversion.cpp version.f90 main.f90

kestrel_FCFLAGS = -Wall -ffree-line-length-0 -fno-range-check -cpp $(NETCDF4_FFLAGS) $(AM_FCFLAGS) $(DEFS)
kestrel_CXXFLAGS = -Wall $(GDAL_CFLAGS) $(PROJ_CFLAGS) $(AM_CXXFLAGS)
kestrel_LDFLAGS = $(GDAL_DEP_LDFLAGS) $(GDAL_LDFLAGS) $(NETCDF4_LDFLAGS) $(PROJ_LDFLAGS)
kestrel_LDADD = $(FCLIBS) $(GDAL_LDFLAGS) $(NETCDF4_FLIBS) $(PROJ_LDFLAGS) $(LIBS)

if WANT_BOOST
kestrel_LDFLAGS += $(BOOST_LDFLAGS)
kestrel_LDADD += $(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB)
else
kestrel_LDADD += -lstdc++fs
endif
if PARAL
kestrel_LDFLAGS += -fopenmp
endif

if HAVE_JULIA
check:
	cd ../tests && $(JULIA) runall.jl
else
check:
	echo "Error: Working Julia executable is required to run test suite" && exit 1
endif

mostlyclean-local:
	rm -f *.mod core
