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
  
Running simulations
-------------------

Kestrel is run on the command line by specifying a text file containing the
required settings and parameters for the simulations. We have provided a couple
of correctly formatted input files to get you started in the `examples/`
directory.

1D example
^^^^^^^^^^

To begin with, we shall set up a simple simulation of an initial release of
turbulent water, propagating down an inerodible constant slope. Take a look at
the file `examples/Input1d_cap_constslope.txt`. You should see some settings
fields divided into conceptual blocks. Let's go through these in turn.

First is the ``Domain`` block, which sets up the simulation area. 

.. literalinclude:: ../../examples/Input1d_cap_constslope.txt
   :lines: 1-8

The domain is a rectangular region centred at the coordinates specified by Lat
and Lon. It is split divided into 6 tiles in the :math:`x` direction. Each of these
tiles is further subdivided into 200 finite volumes. Since this is a 1D
simulation, the corresponding parameters in the :math:`y` direction are set to 1. The
final variable sets the width in metres of each tile. This implies that the
resolution of the numerical discretisation is :math:`\Delta x = 200 / 200 = 1`
metres.

Next, the ``Cap`` block sets the initial release of flowing material.

.. literalinclude:: ../../examples/Input1d_cap_constslope.txt
   :lines: 10-16

The ``Topog`` block describes the initial topographic surface over which the flow
propagates.

.. literalinclude:: ../../examples/Input1d_cap_constslope.txt
   :lines: 18-20

In this case, a constant slope with gradient -0.04 (in :math:`x`). ``Topog
params`` is a list that can specify multiple parameters (if required) for more
complicated topographies.

The ``Parameters`` block declares settings that affect the model equations, such
as certain constants, closure functions and other choices.

.. literalinclude:: ../../examples/Input1d_cap_constslope.txt
   :lines: 22-25

Here, we specified that we want to use the commonplace Chézy (or turbulent
fluid) drag law, give the drag coefficient, 0.04, and turn erosion off.

Next are the ``Solver`` settings, which specify further arguments for the
numerical method.

.. literalinclude:: ../../examples/Input1d_cap_constslope.txt
   :lines: 27-29

In particular this includes the length of time to simulate for (100 seconds).
The last line provides the `CFL condition
<https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition>`_
for the solver, which is used to adaptively limit the time step. This normally
defaults to 0.25, but in 1D, a higher value can be used.

Finally, the ``Output`` block requests 4 simulation outputs. These will be
separated by equal time intervals of :math:`\Delta T = 100 / 4` seconds and
placed in the specified directory.

.. literalinclude:: ../../examples/Input1d_cap_constslope.txt
   :lines: 31-33

We can run the simulation from within the `examples/` directory using the
command

.. code-block:: bash

   $ ../src/kestrel Input1d_cap_constslope.txt

You can find the output in `examples/1d_cap_constslope/`.

Viewing the output
^^^^^^^^^^^^^^^^^^

TODO

2D example
^^^^^^^^^^

Below is a much more involved, annotated Kestrel input file, which specifies a
2D simulation to run on an SRTM topographic map:

.. literalinclude:: InputExample_SRTM.txt
   :language: bash

.. note::

    Input files allow many other settings and parameters to be set than shown in InputExample_SRTM.txt.
    For unspecified settings, default values are used.
    See :ref:`settings_and_parameters` for more details.

.. note:: 
    This simulation requires SRTM topographic data, available from `NASA Earth Data <https://www.earthdata.nasa.gov/>`_ and `USGS Earth Explorer <https://earthexplorer.usgs.gov/>`_. This should be stored locally in an SRTM archive.

The input file is passed to Kestrel as a command line argument:

.. code-block:: bash

  $ ./Kestrel InputExample_SRTM.txt.txt


.. _quick_view:

Viewing the output
^^^^^^^^^^^^^^^^^^

Kestrel can output simulation results as headed-column delimited text files (.txt) and/or NetCDF (.nc) files.  The choice is specified in the input file in the *Outputs* block (see :ref:`settings_and_parameters`).

NetCDF files produced by Kestrel are georeferenced when simulations are performed on a topographic map, so can be opened directly in QGIS to view as a map.

The NetCDF files contain physical and computational variables.  To assist with mapping the critical physical variables, Kestrel NetCDF files are compatible with the `Lahar Flow Map Tools <https://bitbucket.org/markwoodhouse/laharflow_maptools/>`_ QGIS plugin.

An example of the output produced by Kestrel by the InputExample_SRTM.txt input file, post-processed in QGIS (using *Lahar Flow Map Tools*) is shown below.

.. image:: OutputExample_SRTM.png
   :width: 100%
   :align: center

