.. _openmp:

OpenMP Parallelization
======================

OpenMP is a widely used, portable framework for shared-memory parallel programming.  In Kestrel, OpenMP is used to accelerate computationally intensive routines by distributing work across available threads, improving performance for suitable workloads.

.. WARNING::
    OpenMP parallelization in Kestrel is experimental and under-construction.
    There may be changes to the settings during development.

OpenMP support in Kestrel is a work-in-progress, and current scaling with increasing core count is limited. If you are running ensembles of simulations, you will achieve overall better performance by running multiple serial jobs concurrently rather than relying solely on OpenMP parallelization.

Dependencies
------------
The OpenMP implementation in Kestrel requires OpenMP version 5+.

.. _openmp_installation:

Installation
------------

Kestrel's experimental OpenMP implementation is available in the "openmp" branch of the `Kestrel GitHub repository <https://github.com/jakelangham/kestrel/>`_. If you have already cloned the repository, you can switch to this branch using

.. code-block:: bash

  $ git fetch --all
  $ git checkout openmp
  $ git pull

then update the build scripts

.. code-block:: bash

  $ autoreconf -fi


OpenMP is an optional feature that is enabled by passing the flag ``--enable-parallel`` to the configure script.  
You can rename the executable using the ``--program-suffix`` option in the configure script (see :ref:`installation` for additional configuration flags). 
Here, it is used to create a distinctive OpenMP-enabled Kestrel executable called ``kestrel_paral``.

.. code-block:: bash

  $ ./configure --enable-parallel --program-suffix="_paral"

To compile, while cleaning previous compilation outputs

.. code-block:: bash

  $ make clean all
  $ make install

If successful, this creates the executable ``bin/kestrel_paral``.

.. _openmp_run:

Running OpenMP-enabled simulations
----------------------------------

Kestrel with OpenMP can be run without changing input files. The OpenMP settings are controlled through environment variables; see the `OpenMP 5.0 Environment Variables Documentation`_ for details.

.. _OpenMP 5.0 Environment Variables Documentation: https://www.openmp.org/spec-html/5.0/openmpch6.html

Two environment variables are of particular importance:

.. role:: bash(code)
    :language: bash

.. role:: csh(code)
    :language: csh

OMP_NUM_THREADS
    This sets the number of threads to use for parallel regions.

    **Example:**
    set the number of threads to 8.

    On Linux, in a Bash shell:

    .. code-block:: bash

       export OMP_NUM_THREADS=8

    In a C shell:

    .. code-block:: csh

       setenv OMP_NUM_THREADS 8

      

OMP_SCHEDULE
    This controls the schedule kind, optional modifier, and chunk size of all OpenMP-parallelized loops in Kestrel, set using the format ``OMP_SCHEDULE=[modifier,]schedule[,chunk_size]``. 
    
    The modifier (such as ``monotonic`` or ``nonmonotonic``) is optional and can be used to further control how iterations are assigned to threads.

    The choice of scheduling *kind* can significantly affect parallelization performance. The available schedule kinds are [#OpenMP-guide-scheduling]_:

    - **static**:  
      iterations are divided into chunks of size ``chunk_size``, and the chunks are assigned to the threads in the team in a round-robin fashion in the order of the thread number.

    - **dynamic**:  
      iterations are distributed to threads in the team in chunks. Each thread executes a chunk of iterations, then requests another chunk, until no chunks remain to be distributed.

    - **guided**:  
      iterations are assigned to threads in the team in chunks. Each thread executes a chunk of iterations, then requests another chunk, until no chunks remain to be assigned.

    - **auto**:  
      the decision regarding scheduling is delegated to the compiler and/or runtime system. The user gives the implementation the freedom to choose any possible mapping of iterations to threads in the team.

    .. tip::
        For long-running simulations, it is recommended to use the ``auto`` schedule, as it allows the runtime to dynamically adapt the scheduling strategy based on workload and system resources, which can lead to improved performance and efficiency over the course of the simulation.

    See the `OpenMP schedule page <https://www.openmp.org/spec-html/5.0/openmpse49.html#x288-20520006.1>`_ for details of the scheduling options.

    **Example:**
    Set the schedule to dynamic with a chunk size of 4.

    On Linux, in a Bash shell:

    .. code-block:: bash

       export OMP_SCHEDULE="dynamic,4"

    In a C shell:

    .. code-block:: csh

       setenv OMP_SCHEDULE "dynamic,4"

    **Example:**
    Set the schedule to auto (recommended for long simulations):

    On Linux, in a Bash shell:

    .. code-block:: bash

       export OMP_SCHEDULE="auto"

    In a C shell:

    .. code-block:: csh

       setenv OMP_SCHEDULE "auto"

Testing
-------

If you want to check that OpenMP enabled Kestrel is working, the Julia test
suite includes an optional subset of tests for the parallel executable.

.. important::
    The parallel tests require that serial and parallel executables are created as ``bin/kestrel_serial`` and ``bin/kestrel_paral``, respectively.

The tests compare the outputs of the serial and parallel executables to ensure they are identical (within a small numerical tolerance).  The tests also report the speed-up factor of the parallel code.

.. note::
    The tests are designed to be short runtime, and therefore in many cases there is little speed-up from parallelization, and in some cases (particularly for one-dimensional simulations) the OpenMP overhead of parallelizing loops results in slower parallel execution.
    
    For longer running simulations, the benefits of parallelization can become more apparent, and you may observe useful reductions in execution time as the workload increases. 
    
    It is recommended to benchmark your simulations and experiment with different OpenMP scheduling options to achieve optimal parallel performance.

To run the parallel tests, change to the ./tests subdirectory and run:

.. code-block:: bash

   $ ./julia runall.jl parallel

Note that this assumes the existence of a valid Julia symlink within ./tests.

.. [#OpenMP-guide-scheduling] modified from `OpenMP Application Programming Interface Version 5.0 November 2018 <https://www.openmp.org/wp-content/uploads/OpenMP-API-Specification-5.0.pdf>`__.