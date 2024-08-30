# Kestrel #

[![DOI](https://joss.theoj.org/papers/10.21105/joss.06079/status.svg)](https://doi.org/10.21105/joss.06079)

Kestrel is a Fortran code for simulating sediment-laden Earth surface flows,
such as debris flows, landslides and flash flooding.

## Features ##

The software is designed to be easy to install and use. Simulations approximate
solutions to a underlying set of shallow-layer equations that describe a generic
flowing mixture of fluid and sediment. A modular design allows a variety of popular 
closures to be selected at runtime to specialise simulations for particular purposes.

* Fully conservative well-balanced positivity-preserving finite volume solver.
* Dynamic erosion and deposition of sediment.
* Eddy viscosity implementation to ensure well-posedness of morphodynamic flows.
* User-settable basal drag and morphodynamics parametrisations.
* Simulation on user-specified topographies: simple geometric surfaces and digital elevation maps (DEMs).
* Variety of boundary conditions for initiating flows, including from multiple locations.
* Output via geo-referenced NetCDF or plain-text.

An overview of the physical modelling framework may be found 
[here](https://kestrel-unibristol.readthedocs.io/en/latest/physical%20model.html).
For full details, including the numerical implementation, consult our 
paper [[1]](https://arxiv.org/abs/2306.16185). The issue of well-posedness referred to above 
is explored in detail in refs [[1]](https://arxiv.org/abs/2007.15989) and 
[[2]](https://arxiv.org/abs/2306.16185).

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
   to make sure you have automake, autoconf (etc) installed and run

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
this installed, it may be run via `make check`, or by  changing to the ./tests
directory and issuing the command

> `julia runall.jl`

This runs a sequential battery of various tests and can take some time. Unless
you are modifying the code, or are very keen to check that it's working as
expected, you probably do not need to bother with this.

## Documentation ##

More detailed user documentation is available
[here](https://kestrel-unibristol.readthedocs.io/en/latest/index.html).

## Contributing ##

Fixes, improvements and suggestions/bug reports are most welcome. The best way to 
contribute is publically, via the 
[issue tracker](https://github.com/jakelangham/kestrel/issues), though you may also
use the details below to contact us directly. A little guidance on the sorts of
contributions that might be valuable is given [here](https://kestrel-unibristol.readthedocs.io/en/latest/contributing.html).

## Contact ##

Kestrel is developed and maintained by Mark J. Woodhouse
(mark.woodhouse@bristol.ac.uk) and Jake Langham (J.Langham@bristol.ac.uk),
University of Bristol.

## References

* [1] Langham J, Woodhouse MJ, Hogg AJ, Jenkins LT, Phillips JC. 2023 Simulating shallow morphodynamic flows on evolving topographies. _arXiv_ **2306.16185**
* [2] Langham J, Woodhouse MJ, Hogg AJ, Phillips JC. 2021 Linear stability of shallow
morphodynamic flows. _J. Fluid Mech._ **916**.
