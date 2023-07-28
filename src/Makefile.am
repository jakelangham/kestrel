if DEBUG
AM_FCFLAGS = -g -fcheck=all -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=all -Wuninitialized -finit-real=snan -DDEBUG_TIMESTEP -DDEBUG_NEVATIVE_DEPTH -DEBUG_NEGATIVE_CONC -DDEBUG_EXCESS_CONC -DDEBUG_SPD
AM_CXXFLAGS = -g
else
AM_FCFLAGS = -O2 -falign-loops=16 -march=native -mtune=native
AM_CXXFLAGS = -O2 -falign-loops=16 
endif

bin_PROGRAMS = kestrel
include_HEADERS = RasterData.h cUTM.h
kestrel_SOURCES = cUTM.cpp RasterData.cpp SetPrecision.f90 \
	Messages.f90 varStringClass.f90 utilities.f90 Interp2d.f90 \
	utm.f90 RunSettings.f90 GeoTiffRead.f90 DomainSettings.f90 \
	InitConds.f90 Closures.f90 TopogFuncs.f90 Parameters.f90 Limiters.f90 \
	SolverSettings.f90 OutputSettings.f90 TopogSettings.f90 Input.f90 Grid.f90 \
    Equations.f90 HydraulicRHS.f90 MorphodynamicRHS.f90 Redistribute.f90 \
    dem.f90 UpdateTiles.f90 SetSources.f90 NetCDFUtils.f90 \
    Output.f90 Restart.f90 TimeStepper.f90 main.f90

kestrel_FCFLAGS = -Wall -ffree-line-length-0 -fno-range-check -cpp $(NETCDF4_FFLAGS) $(AM_FCFLAGS) $(DEFS)
kestrel_CXXFLAGS = -Wall $(GDAL_CFLAGS) $(AM_CXXFLAGS)
kestrel_LDFLAGS = $(NETCDF4_LDFLAGS)
kestrel_LDADD = $(FCLIBS) $(GDAL_LDFLAGS) $(NETCDF4_FLIBS) $(LIBS)

if WANT_BOOST
kestrel_LDFLAGS += $(BOOST_LDFLAGS)
kestrel_LDADD += $(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB)
else
kestrel_LDADD += -lstdc++fs
endif

if HAVE_JULIA
all-local:
	ln -sf $(JULIA) ../tests/julia
endif

mostlyclean-local:
	rm -f *.mod core
