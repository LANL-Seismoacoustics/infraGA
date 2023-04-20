.. _utilities:

=====================================
Utility Functions
=====================================

In addition to the various analysis and visualization capabilities of InfraPy, several utility functions are also included for a variety of tasks.  These include estimating arrival times for source, using fk analysis results to compute the best beam waveform and spectrum, computing the celerity for a signal from a known source, and writing waveform data from an FDSN or database source to local files.  Information on the various utility functions can be summarized using, :code:`infrapy utils --help`, which will print to screen the utility usage information:

    .. code-block:: none

      Usage: infraga utils [OPTIONS] COMMAND [ARGS]...

        infraga utils - utility functions for infraga usage

      Options:
        -h, --help  Show this message and exit.

      Commands:
        build-g2s-grid   Build grid for range dependent analysis from G2S files
        extract-ecmwf    Extract atmospheric information from an ECMWF netCDF file
        extract-terrain  Extract a line or grid of terrain information

**Building a Range Dependent Grid**

The range dependent methods in infraGA/GeoAc require specifying a grid of nodes and atmospheric structure at each defined by a single atmospheric specification file.  Such information is ingested as, :code:`--atmo-prefix profs/example --grid-lats profs/example_lat.dat --grid-lons profs/example_lon.dat` for the example grid included with the software.  Building a custom grid from a set of G2S specification files can be accomplished using the :code:`inferaga utils build-g2s-grid` utility.  

  .. code-block:: none 

    Usage: infraga utils build-g2s-grid [OPTIONS]

      Construct the numbered specifications and grid files needed to run -rngdep
      methods. Assumes file format from https://g2s.ncpa.olemiss.edu/ (e.g.,
      g2stxt_2022011506_-3.0000_-54.0000.dat)

      Inclusion of source info (location and time), an estimated celerity, and
      profiles at multiple reference times enables construction of a temporally
      varying grid where a node at a distance, r, from the source has an estimated
      time delay of, dt = r / cel, and uses the appropriate atmospheric
      information

      Examples:
          infraga utils build-g2s-grid --g2s-path g2s_dir/ --output-path grid/g2s_grid
          infraga utils build-g2s-grid --g2s-path g2s_dir/ --output-path grid/g2s_grid --src-info '[-20.56989, -175.379975, 2022-01-15T04:14:45]' --celerity-est 0.29

    Options:
      --g2s-path TEXT       Path to G2S specifications
      --output-path TEXT    Output dir + label
      --src-info TEXT       Source info (lat, lon, time) (optional)
      --celerity-est FLOAT  Celerity estimate [km/s] (optional)
      -h, --help            Show this message and exit.

The method parses the latitude, longitude, and date-time information from the file names in a specified directory (:code:`--g2s-path`) and identifies the unique values.  It then cycles through the grid in the appropriate order (that which infraGA/GeoAc's range-dependent method expect) and copies atmospheric files into a numbered set of files.  The below example is used during the LSECE analysis (see Blom, 2023) to build a 2 degree resolution grid in the western US (atmospheric data was pulled from the `NCPA G2S server <http://g2s.ncpa.olemiss.edu/>`_)

  .. code-block:: none
    
    infraga utils build-g2s-grid --g2s-path Artemis_g2s/ --output-path test/test

  .. code-block:: none

    Parsing file list to determine grid and available datetimes...
    File Summary:
      Unique latitudes: [32. 34. 36. 38. 40. 42.]
      Unique longitudes: [-117. -115. -113. -111. -109. -107.]
      Unique times: ['2020-10-27T12:00:00']

    Mapping files to grid nodes...
      Writing atmosphere Artemis_g2s/g2stxt_2020102712_32.0000_-117.0000.dat  -->  test/test.0.met
      Writing atmosphere Artemis_g2s/g2stxt_2020102712_32.0000_-115.0000.dat  -->  test/test.1.met
      Writing atmosphere Artemis_g2s/g2stxt_2020102712_32.0000_-113.0000.dat  -->  test/test.2.met
      Writing atmosphere Artemis_g2s/g2stxt_2020102712_32.0000_-111.0000.dat  -->  test/test.3.met
      Writing atmosphere Artemis_g2s/g2stxt_2020102712_32.0000_-109.0000.dat  -->  test/test.4.met
      Writing atmosphere Artemis_g2s/g2stxt_2020102712_32.0000_-107.0000.dat  -->  test/test.5.met
      Writing atmosphere Artemis_g2s/g2stxt_2020102712_34.0000_-117.0000.dat  -->  test/test.6.met
      Writing atmosphere Artemis_g2s/g2stxt_2020102712_34.0000_-115.0000.dat  -->  test/test.7.met
      Writing atmosphere Artemis_g2s/g2stxt_2020102712_34.0000_-113.0000.dat  -->  test/test.8.met
      ...
      Writing atmosphere Artemis_g2s/g2stxt_2020102712_42.0000_-109.0000.dat  -->  test/test.34.met
      Writing atmosphere Artemis_g2s/g2stxt_2020102712_42.0000_-107.0000.dat  -->  test/test.35.met

    Finished grid construction.
    Run propagation simulations using:
      infraga sph prop --atmo-prefix test/test. --grid-lats test/test.lats.dat --grid-lons test/test.lons.dat

      
In addition to building a standard grid, the methods allow one to include atmospheric models across a range of date-times and specify a source time and location.  Using the distances from the source to various grid nodes and a reference celerity (default of 0.29 km/s, but accessible as :code:`--celerity-est`) a time-varying grid can be constructed so that the atmospheric structure further from the source use a later reference time for the atmosphere.  Such considerations aren't overly important for regional propagation of 100's of kms, but become notable when considering global scale propagation.

**Extracting Terrain Profiles**

Terrain files for use in propagation simulations require some specified geometry (lines for the 2D methods, Cartesian grids for the 3D methods, and latitude/longitude grids for spherical atmospheric layer methods).  The built-in utility for generating such terrain files downloads the `ETOPO1 model file <https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/ETOPO1_Ice_g_gmt4.grd.gz>`_ and extracts the appropriate geometry.  Usage of the method is summarized below.

  .. code-block:: none 

    Usage: infraga utils extract-terrain [OPTIONS]

      Extract lines or grids of terrain information from an ETOPO1 file

      Examples:
          infraga utils extract-terrain --geom line --lat1 40.0 --lon1 -102.5 --azimuth -90.0 --range 750.0 --output-file line_topo.dat
          infraga utils extract-terrain --geom pnt2pnt --lat1 40.0 --lon1 -102.5 --lat2 40.0 --lon2 -110.0 --output-file line_topo.dat
          infraga utils extract-terrain --geom xy-grid --lat1 35.0 --lon1 -110.0 --lat2 45.0 --lon2 -100.0 --lat-ref 40.0 --lon-ref -105.0 --output-file xy_topo.dat
          infraga utils extract-terrain --geom latlon-grid --lat1 35.0 --lon1 -110.0 --lat2 45.0 --lon2 -100.0 --output-file sph_topo.dat

    Options:
      --geom TEXT              Geometry option ('line', 'pnt2pnt', 'xy-grid' or 'latlon-grid')
      --lat1 FLOAT             Latitude of first point (starting point for 'pnt2pnt', lower-left corner for grids)
      --lon1 FLOAT             Longitude of first point (starting point for 'pnt2pnt', lower-left corner for grids)
      --lat2 FLOAT             Latitude of second point (end point for 'pnt2pnt', upper-right corner for grids)
      --lon2 FLOAT             Longitude of second point (end point for 'pnt2pnt', upper-right corner for grids)
      --ref-lat FLOAT          Reference latitude of second point (0.0 for xy-grid option)
      --ref-lon FLOAT          Reference longitude of second point (0.0 for xy-grid option)
      --azimuth FLOAT          Azimuth of great circle path for line option
      --range FLOAT            Great circle distance for line option
      --output-file TEXT       Output file
      --show-terrain BOOLEAN   Visualize terrain results
      --rcvr-file TEXT         File containing stations to plot on map
      --offline-maps-dir TEXT  Use directory for offline cartopy maps
      -h, --help               Show this message and exit.

The :code:`line` and :code:`pnt2pnt` geometry options produce a 2-dimensional profile of terrain that can be used in the :code:`infraga 2d` methods as well as in the recently developed terrain methods included in the `NCPAprop pape <https://github.com/chetzer-ncpa/ncpaprop-release>`_ parabolic equation methods.  The :code:`line` method requires a latitude and longitude for the source location as well as the propagation azimuth and maximum range.  The :code:`pnt2pnt` method accepts the start and end points of a great circle path and extracts the terrain along that path.  In each case, the terrain information is extracted into a file for use and also visualized to screen for review (example shown below).

  .. image:: _static/_images/2d_terrain.png
      :width: 600px
      :align: center

The 3D and latitude/longitude grids require specifying the lower-left and upper-right bounding points and extract the terrain within the region.  For the 3D Cartesian method, range and bearing are computed from a reference location (defaults to the lower-left corner).  As with the 2D methods, the output file is written and the terrain is printed to screen for review.  In the latitude/longitude method, the terrain is visualized on a `Cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_ map to show reference borders and coastlines as in the below example.  Note that although the visualization extends below sea level to show oceanic regions, the extracted file snaps values below sea level to zero as the infrasonic signals will reflect from the ocean surface.  For cases in which a visualized is needed showing terrain including a source and receiver network, the reference latitude and longitude values are plotted as a source and the :code:`--rcvr-file` can be included similar to the :code:`infraga plot map` methods.

  .. image:: _static/_images/terrain_extraction2.png
      :width: 600px
      :align: center


**Extracting Atmosphere Data from ECMWF**

A tool for extracting G2S-format atmospheric files from an ECMWF netCDF format file has been developed, but not robustly evaluated.  Usage info is summarized below, but ongoing evaluation and de-bugging of the method is needed. 

.. code-block:: none 

  Usage: infraga utils extract-ecmwf [OPTIONS]

    Extract G2S-format atmospheric file(s) from an ECMWF netCDF format file.

    Note: method needs evaluation with current ECMWF ERA5 sample files (might be 
    out of date)

    Examples:
        infraga utils extract-ecmwf --ecmwf-file EN19110100.nc --option single  --lat1 30.0 --lon1 -120.0 --output-path test.met
        infraga utils extract-ecmwf --ecmwf-file EN19110100.nc --option grid  --lat1 30.0 --lon1 -120.0 --lat2 40.0 --lon2 -110.0 --output-path test_grid

  Options:
    --ecmwf-file TEXT       ECMWF netCDF file
    --option TEXT           Extraction option ('single' or 'grid')
    --lat1 FLOAT            Latitude of first point (latitude for 'single', lower-left corner for 'grid')
    --lon1 FLOAT            Longitude of first point (longitude for 'single', lower-left corner for 'grid')
    --lat2 FLOAT            Latitude of second point (not used for 'single', upper-right corner for 'grid')
    --lon2 FLOAT            Longitude of second point (not used for 'single', upper-right corner for grids)
    --sample_skips INTEGER  Frequency of samples in the grid option (defaults to 1 to keep every node)
    --output-path TEXT      Output file
    -h, --help              Show this message and exit.
