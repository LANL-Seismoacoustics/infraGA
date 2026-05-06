.. _installation:

=====================================
Installation
=====================================

-----------------
Operating Systems
-----------------

InfraGA/GeoAc can currently be installed on machines running newer versions of Linux or Apple OS X.  Installation on a Windows system requires C/C++ compilation as well as an Anaconda Python installation.

------------
Dependencies
------------

**FFTW**

The weakly non-linear waveform calculations in infraGA/GeoAc leverage Fourier transform methods in the FFTW library.  Installation can be done via `the FFTW download page <http://fftw.org/download.html>`_ or by using an OS package manager (e.g., '*apt-get install -y fftw-dev*' or similar on Linux systems or '*brew install fftw*' on OS X using homebrew).

**OpenMPI**

Parallelized calculations of ray paths is enabled in infraGA/GeoAc via a separate OpenMPI implementation.  These methods are optional, but can greatly reduce computation times for larger simulations.  Similar to FFTW above, installation can be completed via `the OpenMPI download page <http://open-mpi.org/software/ompi/v4.1>`_ or using an OS package manager (e.g., '*apt-get install -y openmpi-bin*' or similar on Linux systems or '*brew install open-mpi*' on OS X).

**Anaconda**

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

At this point if you look in the 'infraga/bin/' directory there should a number of binary files (e.g., :code:`infraga-2d`, :code:`infraga-3d`, :code:`infraga-accel-3d`,...)  These are the C/C++ and OpenMPI binaries to run the various simulations and can be called directly if you're familiar with the pre-Anaconda usage of infraGA/GeoAc.

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


-------------------------------------------
Python Geophysics Suite (PyGS) Installation
-------------------------------------------

Infrasound software tools developed by LANL SMEs have become increasing coupled in usage so that having them in a common Python environment is useful.  An in-development Python Geophysics Suite (PyGS) YML file is included in the InfraPy repository that will build an environment and install InfraPy, infraGA/GeoAc, and stochprop from GitHub.  It can be run using the same syntax as above,

.. code-block:: bash

    >> conda env create -f pygs_env.yml

All dependencies will be installed and the LANL Python libraries pulled from GitHub to complete the environment.  To finish setting up, activate the environment and compile the infraGA/GeoAc software,

.. code-block:: bash

    >> conda activate pygs
    >> infraga compile 

---------------------------------------------------------
Python Geophysics Suite (PyGS) Installation - Dev Version
---------------------------------------------------------

Because the PyGS YML file installs via GitHub cloning, it doesn't copy the examples/ directories from the various libraries for demonstration and also doesn't leave the source code easily accessible for any de-bugging or customization.  A separate developer version is also included that requires a few more steps.  Build an instance of the environment with just InfraPy included using the included YML file,

.. code-block:: bash

    >> conda env create -f pygs-dev_env.yml

Next, clone the other repositories if you don't have them,

.. code-block:: bash

    >> git clone https://github.com/LANL-Seismoacoustics/infrapy.git
    >> git clone https://github.com/LANL-Seismoacoustics/stochprop.git

If you have SSH keys set up for GitHub, you can alternately clone as,

.. code-block:: bash
	
    >> git clone git@github.com:LANL-Seismoacoustics/infrapy.git
    >> git clone git@github.com:LANL-Seismoacoustics/stochprop.git

Once the PyGS development environment is built, activate it using :code:`conda activate pygs_dev` and then use pip with the :code:`-e` flag to install InfraPy and stochprop and compile the infraGA/GeoAc ray tracing methods,

.. code-block:: bash

    >> cd /path/to/infrapy
    >> pip install -e .

    >> infraga compile 

    >> cd /path/to/stochprop
    >> pip install -e .

This installation will clone the example directories with all relevant data and also allow you to interact with other :code:`git` branches for customization.

**Testing**

Once the installation is complete, you can test that the InfraGA/GeoAc methods are set up and accessible by first activating the *infraga* environment with:

.. code-block:: bash

    conda activate infraga_env

The infraGA/GeoAc methods have usage summarizes that can be displayed via the :code:`--help` option.  On the command line, run:

.. code-block:: bash

    infraga --help

The usage information should be displayed:

.. code-block:: none

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

Each of the individual methods have usage information (e.g., :code:`infraga sph --help`) that will be discussed in the :ref:`quickstart`.

**Updates**

It's highly recommended that users keep up with ongoing R&D and related bug fixes and updates of the infraGA/GeoAc methods through the `LANL Seismoacoustics github page <https://github.com/LANL-Seismoacoustics/infraGA>`_.  InfraGA/GeoAc is a research code and bug fixes, new features, and improved utility and plotting methods are frequently pushed to the git repository.  It's possible to download a static copy of the software, but updates can more easily be applied through a git instance where the current version can be readily accessed using a :code:`git pull` command.  The dependency installs above (e.g., FFTW, OpenMPI) don't need to be re-installed and the '-e' option used during the infraga environment build keeps any local updates to the Python methods within the environment; therefore, when updates are applied a simple :code:`make clean`, :code:`make`, :code:`make accel` is all that's needed to update the binaries for infraGA/GeoAc.
