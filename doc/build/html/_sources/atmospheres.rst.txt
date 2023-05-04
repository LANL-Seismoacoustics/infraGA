.. _atmospheres:

=====================================
Atmospheric Data
=====================================

InfraGA/GeoAc requires atmospheric information to specify the acoustic sound speed, ambient winds, and atmospheric density in order to compute propagation information.  A climatological model ('examples/ToyAtmo.met') is included with the git repo for some basic simulations and tests as well as a sample Ground-to-Space (G2S) atmospheric specification.  The climatological model is based on the polynomial fit to the US Standard atmosphere with Gaussian wind structures in the lower- and middle atmosphere as well as idealized atmospheric tidal winds in the meso- and thermosphere (see Lingevitch et al., 1999 and Blom & Waxler, 2021).

The format of atmospheric data can be defined using :code:`prof_format` when using C/C++ binaries or Python configuration files and with :code:`--prof-format` in the Python command line interface.  The default atmospheric file format contains columns describing altitude (km), temperature (Kelvin), zonal (east/west) winds (m/s), meridional (north/south) winds (m/s), density (g/cm :math:`^3`), and pressure (mbar).  The zonal and meridional winds are typically denoted as 'u' and 'v', so that this format is defined as :code:`zTuvdp`.  From this ingested atmospheric data, acoustic sound speed is computed using the adiabatic assumption, :math:`p \propto \rho^\gamma \, \rightarrow \, c^2 = \frac{\partial p}{\partial \rho} = \gamma \frac{p}{\rho}` where :math:`p` is pressure, :math:`\rho` is density, and :math:`\gamma` is the specific heat ratio for the atmosphere (1.4 for dry air).

Two other profile formats are accepted in the current implementation of infraGA/GeoAc: :code:`zuvwTdp` and :code:`zcuvd`.  The first of these also includes a vertical wind component and was used in earlier output options of G2S.  The later defines the sound speed directly and doesn't include temperature or pressure.  The units on these alternate formats are the same as the :code:`zTuvdp` and summarized below in the table for easy reference.

+--------------------------+-------------------------------------------------------------------------------+
| **Profile Format**       | **Fields/Units**                                                              |
+--------------------------+-------------------------------------------------------------------------------+
| :code:`zTuvdp` (default) | z (km), T (K), u (m/s), v (m/s), d :math:`(\text{g/cm}^3)`, p (mbar)          |  
+--------------------------+-------------------------------------------------------------------------------+
| :code:`zcuvd`            | z (km), c (m/s), u (m/s), v (m/s), d :math:`(\text{g/cm}^3)`                  |   
+--------------------------+-------------------------------------------------------------------------------+
| :code:`zuvwTdp`          | z (km), u (m/s), v (m/s), w (m/s), T (K), d :math:`(\text{g/cm}^3)`, p (mbar) |
+--------------------------+-------------------------------------------------------------------------------+

********************
Climatological Model
********************

A climatological model is included for some general analysis and leveraged for a number of the :ref:`quickstart` examples.  This file is located in 'examples/ToyAtmo.met' and includes the polynomial fit introduced by Lingevitch et al. (1999) with Gaussian wind structures corresponding to an eastward jet stream near the tropopause, westward stratospheric wind jet (typical of summer atmospheric structure in the northern hemisphere), and sinusoidal atmospheric tides.  Detailed specifications for these winds are summarized in Blom & Waxler (2021).  The below figure shows the atmospheric structure and waveguide refraction heights using the :code:`infraga plot atmo` functionality discussed as discussed during the :ref:`quickstart`.

  .. image:: _static/_images/atmo_plot.png
      :width: 1000px
      :align: center

**************************
Atmospheric Data Resources
**************************

The Ground-to-Space (G2S) methodologies (Drob et al., 2003) interpolate and merge numerical weather prediction data with climatological models for the meso- and thermosphere to produce atmospheric structure information that extends from the ground surface to above 100 km altitude.  Currently, a server hosting G2S atmospheric data is run by infrasound experts at the  
`University of Mississippi's National Center for Physical Acoustics <https://g2s.ncpa.olemiss.edu/>`_.  Through this web interface, one can request atmospheric data files for single locations, along great-circle paths, or on latitude/longitude grids.  The interface for the server is relatively straightforward to navigate and requested atmospheric data is shared via a download link sent to the requester's email.

  .. image:: _static/_images/G2S_NCPA.png
      :width: 900px
      :align: center


Atmospheric data files from the NCPA G2S system have file names that specify the date and time of the atmospheric information as well as its location on the globe (e.g., `'g2stxt_2020102922_32.0000_-107.0000.dat'`).  Each file includes a header of information summarizing the model data source, when it was constructed, the reference time and location as well as the list of fields included and the ground elevation at the location.  Also included are a series of line formatted to be read by the `NCPAprop software <https://github.com/chetzer-ncpa/ncpaprop-release>`_.  This header information is not utilized by infraGA/GeoAc and the file format must be specified as part of the simulation as noted above (note: the current G2S file format is the default infraGA/GeoAc format, so in general no modifications are necessary unless using another data source).  An example NCPA G2S file header is below.

  .. code:: none

    # Data Source: NASA MERRA (version 2)
    # Model Calculated 2020-11-28 22:44:02
    # Model Time = 2020-10-29 22:00:00
    # Location = [  32.0000, -107.0000 ]
    # Fields = [ Z(km), T(K), U(m/s), V(m/s), R(g/cm3), P(mbar) ]
    # Ground Height = 1.292
    # The following lines are formatted input for ncpaprop
    #% 0, Z0, km, 1.292
    #% 1, Z, km
    #% 2, T, degK
    #% 3, U, m/s
    #% 4, V, m/s
    #% 5, RHO, g/cm3
    #% 6, P, mbar
      0.000      2.97874e+02      -3.99578e+00       1.72398e+00       1.19409e-03       1.02082e+03
      0.100      2.97214e+02      -3.97783e+00       1.72879e+00       1.18278e-03       1.00891e+03
      0.200      2.96561e+02      -3.97106e+00       1.75268e+00       1.17155e-03       9.97141e+02
      ...


In addition to the NCPA G2S server, the framework developed by Drob et al. (2003) for G2S construction has been reproduced by infrasound and atmospheric scientists at the Alaska Volcano Observatory (AVO) and made available as `AVO G2S <https://github.com/DOI-USGS/volcano-avog2s>`_.  The methods included require some additional downloads of USGS libraries but provide a more robust system able to ingest additional low- and middle atmosphere weather data; however, for simplicity of use, the NCPA G2S server is the recommended means of requesting and obtaining atmospheric data for infrasonic propagation modeling.

  .. image:: _static/_images/AVO_G2S.png
      :width: 900px
      :align: center

Lastly, the European Centre for Medium-Range Weather Forecasts (ECMWF) provides `ERA5 Reanalysis atmospheric specifications <https://registry.opendata.aws/ecmwf-era5/>`_ in netCDF4 format files.  Preliminary methods to extract G2S-style files for use in infraGA/GeoAc as well as the NCPAprop methods are included in infraGA/GeoAc's :ref:`utilities`.