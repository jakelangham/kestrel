# Kestrel #

Kestrel is a Fortran code for simulating sediment-laden Earth surface flows,
such as debris flows, landslides and flash flooding.

## 'Very' quick start ##

Most likely, you will need to refer to the
[documentation](https://kestrel-unibristol.readthedocs.io/en/latest/index.html)
to use Kestrel effectively and in an informed manner. This includes a proper
[quick start](https://kestrel-unibristol.readthedocs.io/en/latest/quickstart.html)
guide aimed towards new users.

For the particularly impatient, here is the basic information needed to get
started

1. Kestrel builds with GNU autotools on Unix-like platforms. The auto-generated
   configuration scripts are not included in the git repository, so you'll need
   to do

   > `autoreconf -fi`

   followed by

   > `./configure && make`

   which places the Kestrel executable in the src/ directory. You will need an
   up-to-date [GCC](https://gcc.gnu.org/) (9+), the [GDAL](https://gdal.org/),
   [PROJ](https://proj.org) and (optionally) the
   [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) libraries.

3. Simulations are conducted by specifying an input file on the command line,
   i.e.

   > `/path-to-kestrel [path-to-input-file]`

   There are some examples of correctly formatted input files included in this
   repository in the examples/ directory. These may be adapted to suit different
   needs.

4. To conduct simulations on measured topographies, for example via digital
   elevation models (DEMs), you will need provide this data and tell Kestrel
   where to look for the appropriate files. See either the example input files,
   or the 
   [user documentation](https://kestrel-unibristol.readthedocs.io/en/latest/index.html)
   for how to do this.

## Testing ##

The full test suite for Kestrel makes use of the Julia language. If you have
this installed, it may be run by changing to the ./tests directory and issuing
the command

> `julia runall.jl`

This runs a sequential battery of various tests and can take some time. Unless
you are modifying the code, or are very keen to check that it's working as
expected, you probably do not need to bother with this.

## Documentation ##

More detailed user documentation is available
[here](https://kestrel-unibristol.readthedocs.io/en/latest/index.html).

## Contact ##

Kestrel is developed and maintained by Mark J. Woodhouse
(mark.woodhouse@bristol.ac.uk) and Jake Langham (J.Langham@bristol.ac.uk),
University of Bristol.
