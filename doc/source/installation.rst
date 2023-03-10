.. _installation:

=====================================
Installation
=====================================

-------------------------------------
Operating Systems
-------------------------------------

InfraGA/GeoAc can currently be installed on machines running newer versions of Linux or Apple OSX.  Installation on Windows system requires C/C++ compilation as well as an Anaconda Python installation.  

-------------------------------------
FFTW
-------------------------------------

The weakly non-linear waveform calculations in infraGA/GeoAc leverage Fourier transform methods in the FFTW library.  Installation can be done via `download <http://fftw.org/download.html>`_ or by using an OS package manager (e.g., '*apt-get install -y fftw-dev*' or similar on Linux systems or '*brew install fftw*' on OS X using homebrew).

-------------------------------------
OpenMPI
-------------------------------------

Parallelized calculations of ray paths is enabled in infraGA/GeoAc via a separate OpenMPI implementation.  These methods are optional, but can greatly reduce computation times for larger simulations.  Similar to FFTW above, installation can be completed via `download <http://open-mpi.org/software/ompi/v4.1>`_ or using an OS package manager (e.g., '*apt-get install -y openmpi-bin*' or similar on Linux systems or '*brew install open-mpi*' on OS X).

-------------------------------------
Anaconda
-------------------------------------

infraGA/GeoAc methods can be called directly from the command line without the Python wrappers; though, most users will find the wrappers, utilities, and visualization tools in the *infraga* Python interface to be easier to use.  Installation of these Python tools currently depends on Anaconda to resolve and download the correct python libraries.  Anaconda cannot currently be installed via package managers and must be `downloaded and installed <https://www.anaconda.com/distribution/>`_ manually.

-------------------------------------
infraGA/GeoAc Installation
-------------------------------------

Installation of infraGA/GeoAc requires several steps to 1) build the normal C/C++ single-CPU version of the software, 2) build the OpenMPI multi-threaded version, and 3) build the Anaconda environment that contains the wrappers, utilities, and plotting functions.

From the primary directory for infraGA/GeoAc (containing several directories: 'bin', 'doc', 'examples', etc. as well as a makefile), simply run

.. code-block:: bash 

    make 
    
and the C/C++ single-CPU version of the software will compile.  Similarly, the accelerated multi-threaded methods can be compiled using,

.. code-block:: bash 

    make accel

At this point if you look in the :code:`infraga/bin/` directory there should a number of binary files (e.g., :code:`infraga-2d`, :code:`infraga-3d`, :code:`infraga-accel-3d`,...)  These are the C/C++ and OpenMPI binaries to run the various simulations and can be called directly if you're familiar with the pre-Anaconda usage of infraGA/GeoAc.

Once the binaries have been compiled, the *infraga* Anaconda environment can be built by running:

.. code-block:: bash

    conda env create -f infraga_env.yml

If this command executes correctly and finishes without errors, it should print out instructions on how to activate and deactivate the new environment:

To activate the environment, use:

.. code-block:: bash

    conda activate infraga_env

To deactivate an active environment, use

.. code-block:: bash

    conda deactivate

Note: in some cases Linux systems require activation via :code:`source activate infraga_env`, so if you're getting an Anaconda error about activating the environment try this alternate method.

-------------------------------------
Testing
-------------------------------------

Once the installation is complete, you can test that the InfraGA/GeoAc methods are set up and accessible by first activating the environment with:

.. code-block:: bash

    >> conda activate infraga_env

The infraGA/GeoAc methods have usage summarizes that can be displayed via the :code:`--help` option.  On the command line, run:

.. code-block:: bash

    infraga --help

The usage information should be displayed:

.. code-block:: bash

    Usage: infraga [OPTIONS] COMMAND [ARGS]...

      infraga - Python interface for using the infraGA/GeoAc software tools

    Options:
      -h, --help  Show this message and exit.

    Commands:
      2d     Run 2d (effective sound speed) ray tracing
      3d     Run 3d (moving medium) ray tracing
      plot   Various visualization functions
      sph    Run spherical layer (moving medium) ray tracing
      utils  Various utility functions

Each of the individual methods have usage information (e.g., :code:`infraga 3d --help`) that will be discussed in the :ref:`quickstart`

