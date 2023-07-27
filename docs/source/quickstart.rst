.. _quick_start:

Quick start
===========

.. _dependencies:

Dependencies
------------
Kestrel is written in Fortran 2003, with some additional C++ calls to geospatial libraries.
It should be possible to compile and run it on most modern Unix-like platforms relatively painlessly.
In order to do so, make sure the following dependencies are present on your system:

* Either `GCC <https://www.gnu.org/software/gcc/>`_ version 9+, or an older GCC plus the `Boost <https://www.boost.org/>`_ C++ libraries (>= 1.46.0)
* `GDAL <https://gdal.org/>`_ (>= 2.2.0)
* `PROJ <https://proj.org/>`_
* GNU autotools (if building directly from the git repository)
* (optional) `NetCDF <https://www.unidata.ucar.edu/software/netcdf/>`_
* (optional) `Julia <https://julialang.org/>`_, for running tests

.. _installation:

Installation
------------

We use GNU autools for our build system. If building directly from the git repository, you must first generate the relevant build scripts using the command

.. code-block:: bash

  $ autoreconf -fi

Then the code can be compiled with

.. code-block:: bash

  $ ./configure && make

If succesful, this places an executable in the src directory.

You may wish to provide additional flags to the configure script, depending on
your local setup. For example, you can do

.. code-block:: bash

  $ ./configure --with-netcdf4=yes/no/PATH

to specify whether to build Kestrel with NetCDF support, where PATH is the
base directory of your (serial) NetCDF installation.  Since NetCDF is useful
(and preferred) we build this support in, if found, unless otherwise specified.

Likewise, if you do not have a newer (>= 9) version of GCC,
you can instead use Boost to do filesystem calls by specifying

.. code-block:: bash

  $ ./configure --with-boost[=ARG]

where ARG optionally specifies the location of your Boost installation.

For testing, the path to a valid Julia executable can be specified with

.. code-block:: bash

  $ ./configure JULIA=[PATH]

This will place a symlink to Julia in the ./tests directory.

Some other options are listed in the help dialogue of the configure script

.. code-block:: bash

  $ ./configure --help

.. _quick_run:
  
Running a simulation
--------------------

Kestrel simulations require a domain, topography, sources of material, specification of model parameters (including choices of closures), output formats, and solver settings. These are passed to Kestrel in a Kestrel Input file, with blocks corresponding to each of the settings.  Details of each of these requirements are provided in :ref:`settings_and_parameters`.

Below is an example, annotated Kestrel Input file

.. literalinclude:: InputExample_SRTM.txt
   :language: bash
   :linenos:

The input file is passed to Kestrel as a command line argument:

.. code-block:: bash

  $ ./Kestrel InputExample.txt


.. _quick_view:

Viewing the output
------------------

Kestrel can output simulation results as headed-column delimited text files (.txt) and/or NetCDF (.nc) files.  The choice is specified in the input file in the *Outputs* block (see :ref:`settings_and_parameters`).

NetCDF files are georeferenced when simulations are 