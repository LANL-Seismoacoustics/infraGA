.. _advanced:

=====================================
Advanced Usage
=====================================

In addition to the :code:`prop` option that simulate propagation from a point source through a specified atmospheric structure.  Additional parameter settings are available to include range dependence and realistic terrain interactions and additional simulation options include methods to compute eigenray paths connecting specific source-receiver geometries as well as compute the weakly non-linear (weak shock limit) waveform evolution along individual ray paths.

****************************
Advanced Parameter Settings
****************************

**Range Dependent Analysis**

Using the :code:`--atmo-file` parameter defines a single atmospheric specification to use for a stratified atmospheric model.  In cases where propagation extends more than a few hundred kilometers, horizontal variations in the atmospheric structure can become significant.  Range dependent atmospheric structure can be used in infraGA/GeoAc via specification of an atmospheric file prefix as well as horizontal grid nodes.  An example 4x4 grid is provided in the examples/profs directory that includes 16 atmospheric files (example0.met, example1.met, example2.met, etc.) along with example_lat.dat and example_lon.dat that contain the latitude and longitude of the grid nodes.  

    .. code:: none

      infraga sph prop --atmo-prefix profs/example --grid-lats profs/example_lat.dat --grid-lons profs/example_lon.dat --src-lat 40.0 --src-lon -102.5 --azimuth -90.0 --z-grnd 1.0 


    .. code:: none

      	############################################
        ####     Running infraga-sph-rngdep     ####
        ####             Propagation            ####
        ############################################

      Interpolating atmosphere data in 'profs/example'* using format 'zTuvdp'...
        Setting grid node at (35, -110) with profile profs/example0.met
        Setting grid node at (35, -106.67) with profile profs/example1.met
        Setting grid node at (35, -103.33) with profile profs/example2.met
        ...
        Setting grid node at (38.33, -110) with profile profs/example4.met
        Setting grid node at (38.33, -106.67) with profile profs/example5.met
        Setting grid node at (38.33, -103.33) with profile profs/example6.met
        ...
        Setting grid node at (41.66, -110) with profile profs/example8.met
        Setting grid node at (41.66, -106.67) with profile profs/example9.met
        Setting grid node at (41.66, -103.33) with profile profs/example10.met
        ...

        Propagation region limits:
          latitude = 35.001, 44.999
          longitude = -109.999, -100.001
          altitude = 0, 150


      Parameter summary:
        inclination: 0.5, 45, 0.5
        azimuth: -90, -90, 1
        bounces: 2
        source location (lat, lon, alt): 40, -102.5, 1
        ground elevation: 1
        frequency: 0.1
        S&B atten coeff: 1
        write_atmo: false
        write_rays: true
        write_topo: false
        write_caustics: false
        calc_amp: true

      Calculating ray path: 0.5 degrees inclination, 90 degrees azimuth.
      Calculating ray path: 1 degrees inclination, 90 degrees azimuth.
      ...
      Calculating ray path: 44.5 degrees inclination, 90 degrees azimuth.
      Calculating ray path: 45 degrees inclination, 90 degrees azimuth.

Note that the various atmospheres are read in and set on the grid nodes in a specific order (cycling through longitude values first).  When using a user created grid it's useful to check that the atmospheric specifications are being ingested and set at the correct nodes.  A utility function is available to build a grid for range dependent analysis in :ref:`utilities`.

Visualization of range-dependent results is relatively straightforward and similar those doing so with stratified results; though, only one of the atmospheric files can be visualized so that the arrivals and ray path files must be specified directly.

  .. code:: none
    
    infraga plot azimuthal --atmo-file profs/example0.met --arrivals profs/example.arrivals.dat --ray-paths profs/example.raypaths.dat

  .. image:: _static/_images/rng-dep_azimuthal.png
      :width: 1200px
      :align: center


Although the range-dependent effects aren't overly obvious in this simulation due to the limited propagation range, comparison of the inclination angles and point-to-point propagation of the ray paths does exhibit variations.  Due to the added complexity of multi-variate interpolation needed to compute the gradients of the atmospheric parameters not only in the vertical direction but also in the horizontal, computation of range-dependent ray paths is notably slower than stratified analyses.

For global scale analyses (e.g., analysis of propagation paths from large bolides or volcanic eruptions), one can define a maximum propagation range greater than half of the circumference of the Earth (~24,000 km) so that the great circle distances never reach the break condition and use the :code:`--max-tm` parameter to define a maximum propagation time (defined in hours): :code:`--max-rng 25000 --max-tm 24.0 --bounces 1000`

**Including Realistic Terrain**

The simulations run in the :ref:`quickstart` assumed a flat ground that could be modified using :code:`--z-grnd` to change its elevation; however, infrasonic signals have wavelength that can be comparable in scale to large scale terrain features (e.g., mountain ranges, canyons) so that interaction of energy with realistic terrain is an important factor in propagation simulations.  A terrain file can be specified for use in a simulation using the parameter :code:`--topo-file` (or :code:`topo_file` in the C/C++ or config file syntax).  For a spherical geometry simulation, this file should contain the elevation of the ground surface (in kilmoeters relative to sea level) at a series of latitude and longitude points.  One of the included :ref:`utilities` can be used to build a grid file for use:

  .. code-block:: none

    infraga utils extract-terrain --geom latlon-grid --lat1 35.0 --lon1 -110.0 --lat2 45.0 --lon2 -100.0 --output-file sph_topo.dat

In this function call a latitude/longitude grid is defined with a lower-left corner at (35.0, -110.0) and upper-right corner at (45.0, -100.0).  The terrain is written into a file named :code:`sph_topo.dat` and a map showing the terrain is printed to screen.  The method will download the `ETOPO1 model file <https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/ETOPO1_Ice_g_gmt4.grd.gz>`_ for use and store it within the infraGA/GeoAc directory structure at an expected location for repeated use in subsequent function calls.  The spatial resolution of the output preserves the 1 arc-minute (:math:`\sim1.85` km) resolution of the ETOPO1 model.  Also, it should be noted that the elevation scale extends below sea level to show oceanic regions in the visualization; however, the terrain file snaps such locations to sea level as infrasonic waves will reflect from the ocean surface.

  .. image:: _static/_images/terrain_extraction.png
      :width: 600px
      :align: center


The file created by the utility function contains 3 columns describing latitude, longitude, and ground surface elevation.  The ordering of grid points is such that longitude values are cycled through at each latitude point (note that this is the same convention used in specifying files for the range-dependent grid above).

  .. code-block:: none

    35.00000000000003 -110.0              1.639
    35.00000000000003 -109.98333333333333 1.62
    35.00000000000003 -109.96666666666667 1.634
    ...
    35.00000000000003 -100.01666666666667 0.607
    35.00000000000003 -100.0              0.582
    35.016666666666694 -110.0             1.657
    35.016666666666694 -109.98333333333333 1.648
    35.016666666666694 -109.96666666666667 1.649
    ...

The terrain file can used in a simulation by specifying it as the topography file (:code:`--topo-file` or :code:`topo_file`).  It should be noted that given the limited latitude and longitude bounds of the grid, the default source location in the :code:`infraga sph` methods may not be valid and the default source location when user input is not provided is defined by the mid-point of the propagation domain.  In the following simulation, a source is defined to the east of the Rocky Mountains in the western US and ray paths are computed for propagation in the westward stratospheric waveguide:

  .. code-block:: none

    infraga sph prop --atmo-file ToyAtmo.met --topo-file sph_topo.dat --src-lat 40.0 --src-lon -102.5 --azimuth -90.0


  .. code-block:: none

      ###########################################
      ####        Running infraga-sph        ####
      ####            Propagation            ####
      ###########################################

    Interpolating atmosphere data in 'ToyAtmo.met' and topography data in 'sph_topo.dat'...
      Propagation region limits:
        latitude = 35.001, 44.9823
        longitude = -109.999, -100.001
        altitude = 0, 139.9

      Maximum topography height: 4.106
      Boundary layer height: 6.106

    Parameter summary:
      inclination: 0.5, 45, 0.5
      azimuth: -90, -90, 1
      bounces: 10
      source location (lat, lon, alt): 40, -102.5, 1.208
      turn_ht_min: 0.2
      frequency: 0.1
      S&B atten coeff: 1
      write_atmo: false
      write_rays: true
      write_topo: true
      write_caustics: false
      calc_amp: true

    Calculating ray path: 0.544095 degrees inclination, -90 degrees azimuth.
    Calculating ray path: 1.0441 degrees inclination, -90 degrees azimuth.
    Calculating ray path: 1.5441 degrees inclination, -90 degrees azimuth.
    ...

The output is similar to the flat ground simulation; though, the terrain file ingestion is noted and the maximum terrain height and an approximated boundary layer height (2 km above the peak of terrain) is printed to screen for reference.  

The ground intercept condition is more numerically intensive when considering terrain because instead of comparing the ray altitude to a constant value the interpolated ground surface must be evaluated at the latitude and longitude of the ray path.  The reason that a 2 km layer is chosen is that the envelope which smoothly forces the wind fields to zero at the ground surface extends some distance above the ground and that envelope also requires evaluation of the terrain interpolation.  This results in longer computation times for tropospheric propagation paths for which the ray path is in this near ground region for a large portion of the simulation, but avoid increasing the computation time for middle- and upper atmospheric portions of the path when terrain isn't directly impacting the propagation path.

Visualization of the propagation results can be done using the same azimuthal function from the :ref:`quickstart`, :code:`infraga plot azimuthal --atmo-file ToyAtmo.met`, and the interaction of terrain is evident in the scatter of reflected arrivals near 350 km downrange that becomes even more severe beyond 500 km.

  .. image:: _static/_images/terrain_prop1.png
      :width: 1200px
      :align: center

The interaction with terrain can be directly seen by zooming the figure into the near ground as shown below (simply click/drag to highlight the region in the matplotlib window).  The variable altitude and gradients of the terrain are evident in the lower bound of the ray paths.

  .. image:: _static/_images/terrain_prop2.png
      :width: 1200px
      :align: center

It would be useful to also visualize the terrain along the propagation path; however, this isn't immediately possible as the terrain file defined the latitude and longitude grid for the entire region.  The parameter :code:`--write-topo` can be used to output the terrain profile immediately underneath the first ray path.  Though not included in the below example, it is also useful to leverage the OpenMPI methods for this computation because of the increased computation time needed to use the terrain interpolation.  Recalling that the multi-threaded methods automatically turn off ray path output, this requires :code:`--write-rays True --cpu-cnt 8` with the CPU count set to whatever is available on you machine.

  .. code-block:: none

    infraga sph prop --atmo-file ToyAtmo.met --topo-file sph_topo.dat  --azimuth -90.0 --src-lat 40.0 --src-lon -102.5 --write-topo True

This simulation produces the familiar raypaths and arrivals files as well as a file named :code:`ToyAtmo.terrain.dat` that includes the latitude, longitude, and terrain elevation below the first ray path computed.  This terrain profile can be included in the azimuthal visualization to show the terrain along the propagation azimuth,

  .. code-block:: none
    
    infraga plot azimuthal --atmo-file ToyAtmo.met --terrain-profile ToyAtmo.terrain.dat

  .. image:: _static/_images/terrain_prop3.png
      :width: 1200px
      :align: center

This visualization is notably more clear in showing what the terrain profile looks like along the path.  

It should be noted that in some cases cross winds and terrain interactions can lead to ray paths that deviate significantly from the initial azimuthal plane.  In such a case, the terrain profile projected below the first ray path may differ significantly from that below other ray paths and those ray paths may appear to reflect from regions above or below the projected terrain profile.  In such a case, a secondary simulation using a higher inclination initial ray path might be useful to investigate the terrain profile below others in the simulation (be sure to save the arrivals and ray path output before doing this as they will be overwritten in the repeated simulations).


**Elevated Sources**

The default inclination angle range used in ray tracing assumes a source at the ground surface and a slightly positive inclination minimum.  In the case of an elevated source, this results in a large portion of ray paths being left uncomputed.  Consider running a simulation with a source at an altitude of 20 kilometers:


  .. code-block:: none
    
    infraga sph prop --atmo-file ToyAtmo.met --src-alt 20.0


The resulting ray paths can be visualized to see how a large portion of energy is missing from the simulation output:


  .. code-block:: none
    
    infraga plot azimuthal --atmo-file ToyAtmo.met

  .. image:: _static/_images/elevated_source1.png
      :width: 900px
      :align: center

The Jacobian computation that is required to solve the Transport equation and calculate the geometric spreading has a singularity for inclination angles of :math:`\pm 90^o`; therefore, a minimum inclination angle of :math:`89^o` should be used as the lower limit in computing ray paths for an elevated source:

  .. code-block:: none
    
    infraga sph prop --atmo-file ToyAtmo.met --src-alt 20.0 --incl-min -89.0

Visualizing this result shows the additional ray paths which are initially downward propagating and ensonify the region below the source.  Also, the initially upward and downward paths produce alternating ensonification in the stratospheric waveguide:

  .. code-block:: none
    
    infraga plot azimuthal --atmo-file ToyAtmo.met

  .. image:: _static/_images/elevated_source2.png
      :width: 900px
      :align: center


Reducing the maximum propagation range to 400 km (adding :code:`--max-rng 400` to the above simulation command) more clearly shows the ensonification below the source and the two sets of ray paths contained within the stratospheric waveguide.  


  .. image:: _static/_images/elevated_source3.png
      :width: 700px
      :align: center


*********************
Advanced Option Usage
*********************

**Eigenray Analysis**

In some scenarios, those specific propagation paths connecting known source and receiver locations are of interest.  Such propagation paths are termed 'eigenrays' and can be difficult to compute when considering propagation paths in 3 dimensions including cross winds.  The auxiliary parameters that are utilized by infraGA/GeoAc for computation of the geometric spreading losses can also be leveraged for computation of launch angle changes that shift an arrival closer to a desired location.  A Levenberg-Marquardt (LM) algorithm has been implemented that uses the auxiliary parameters for such a search (more detail in the :ref:`physics` discussion and in Blom & Waxler (2017)).  The eigenray methods in infraGA/GeoAc are accessed using the :code:`eigenray` flag instead of :code:`prop` and have a number of common parameters.  In addition to specifying a source location, the receiver location is also needed.  The following command runs an eigenray search for a source at (30, -100) to a receiver to the west-south-west at (30.25, -104.25).  

  .. code::  none

    infraga sph eigenray --atmo-file ToyAtmo.met --src-lat 30.0 --src-lon -100.0 --rcvr-lat 30.25 --rcvr-lon -104.25 --bnc-max 1 --verbose True


  .. code:: none

      ###########################################
      ####        Running infraga-sph        ####
      ####          Eigenray Search          ####
      ###########################################

    Interpolating atmosphere data in 'ToyAtmo.met' using format 'zTuvdp'...
      Propagation region limits:
        latitude = -90, 90
        longitude = -180, 180
        altitude = 0, 139.9


    Parameter summary:
      source location (lat, lon, alt): 30, -100, 0
      receiver location (lat, lon, alt): 30.25, -104.25, 0
      inclination range: 0.5, 45
      bounces: 0, 1
      ground elevation: 0
      damping:0.001
      frequency: 0.1
      S&B atten coeff: 1
      verbose: true

    Searching for 0 bounce eigenrays.
      Estimating eigenray angles for source-receiver at great circle distance 409.604 km and azimuth -85.0442 degrees from N.  Inclination limits: [0.5, 45].
        Ray launched with inclination 0.5, azimuth -85.0442 arrives at range 216.761 km after 0 bounce(s).	Exact arrival at 30.1468 degrees N latitude, -102.247 degrees E longitude
        Ray launched with inclination 0.6, azimuth -85.0442 arrives at range 216.445 km after 0 bounce(s).	Exact arrival at 30.1466 degrees N latitude, -102.243 degrees E longitude
        ...
        Ray launched with inclination 35.8, azimuth -85.0442 arrives at range 410.169 km after 0 bounce(s).	Exact arrival at 30.24 degrees N latitude, -104.256 degrees E longitude
        Ray launched with inclination 35.9, azimuth -85.0442 arrives at range 408.259 km after 0 bounce(s).	Exact arrival at 30.239 degrees N latitude, -104.237 degrees E longitude
        Azimuth deviation = 0.163077, less than 2 degrees: estimates acceptable.

        Searching for exact eigenray using auxiliary parameters.
        Calculating ray path: 35.8 degrees inclination, -85.0442 degrees azimuth		Arrival at (30.240036, -104.25644), distance to receiver = 1.2688308 km.
        Calculating ray path: 35.831669 degrees inclination, -84.874356 degrees azimuth		Arrival at (30.249983, -104.24998), distance to receiver = 0.002455279 km.
        Eigenray-0:
          inclination [deg] = 35.831669
          azimuth [deg] = -84.874356
          bounces [-] = 0
          latitude [deg] = 30.249983
          longitude [deg] = -104.24998
          time [s] = 1503.657
          celerity [km/s] = 0.27240416
          turning height [km] = 130.22471
          arrival inclination [deg] = 35.84037
          back azimuth [deg] = 92.992416
          attenuation (geometric) [dB] = -54.837386
          absorption [dB] = -23.705732

    	Estimating eigenray angles for source-receiver at great circle distance 409.60414 km and azimuth -85.044241 degrees from N.  Inclination limits: [35.9, 45].
        Ray launched with inclination 35.9, azimuth -85.0442 arrives at range 408.259 km after 0 bounce(s).	Exact arrival at 30.239 degrees N latitude, -104.237 degrees E longitude
        ...
        Ray launched with inclination 44.9, azimuth -85.044241 arrives at range 314.57557 km after 0 bounce(s).	Exact arrival at 30.169893 degrees N latitude, -103.26423 degrees E longitude
        Reached maximum inclination angle or iteration limit.

    Searching for 1 bounce eigenrays.
      Estimating eigenray angles for source-receiver at great circle distance 409.60414 km and azimuth -85.044241 degrees from N.  Inclination limits: [0.5, 45].
        Ray launched with inclination 0.5, azimuth -85.044241 arrives at range 433.04978 km after 1 bounce(s).	Exact arrival at 30.255722 degrees N latitude, -104.49409 degrees E longitude
        Ray launched with inclination 0.6, azimuth -85.044241 arrives at range 432.46918 km after 1 bounce(s).	Exact arrival at 30.255473 degrees N latitude, -104.48806 degrees E longitude
        ...

    Identified 3 eigenray(s).

As with the :code:`prop` simulations, atmospheric data is ingested and interpolated to define the propagation medium and the parameter summary provides an overview of the run settings.  

The methodology of infraGA/GeoAc's eigenray search is separated into two stages: In the initial stage, rays are launched in the direction from the source to the receiver at increasing inclination angles.  Once a pair of rays are identified which pass over the receiver range, the LM algorithm is used to search for the exact eigenray.  The search is the resumed from the launch angle that triggered the LM search and these steps are repeated until the maximum inclination angle is reached.  The search then begins with an increased number of ground reflections and continues until the maximum number of such reflections is reached.  Upon completion, the various eigenray paths are written into unique files (e.g., 'ToyAtmo.eigenray-0.dat') and all eigenray arrivals are written into an arrivals file.

In the above simulation, the :code:`--verbose` parameter is set to True and additional information is printed to screen as the eigenray search stages are completed.  If this parameter is set to :code:`false`, then less information is printed to screen as the simulation is run:

  .. code:: none

      ###########################################
      ####        Running infraga-sph        ####
      ####          Eigenray Search          ####
      ###########################################

    Interpolating atmosphere data in 'ToyAtmo.met' using format 'zTuvdp'...
      Propagation region limits:
        latitude = -90, 90
        longitude = -180, 180
        altitude = 0, 139.9


    Parameter summary:
      source location (lat, lon, alt): 30, -100, 0
      receiver location (lat, lon, alt): 30.25, -104.25, 0
      inclination range: 0.5, 45
      bounces: 0, 1
      ground elevation: 0
      damping:0.001
      frequency: 0.1
      S&B atten coeff: 1
      verbose: false

    Searching for 0 bounce eigenrays.
      Eigenray identified:	theta, phi = 35.831669, -84.874356 degrees.
    Searching for 1 bounce eigenrays.
      Eigenray identified:	theta, phi = 4.1392114, -84.970504 degrees.
      Eigenray identified:	theta, phi = 31.703267, -84.583914 degrees.
    Identified 3 eigenray(s).

The verbose option is fairly useful when the user is unsure of what eigenrays are expected or the method is not identify an expected eigenray.  Several other parameters have unique functionality in eigenray searches:

  *Reflection count*: the parameters :code:`--bnc-min` :code:`--bnc-max`, and :code:`--bounces` allow the user to define the range possible ground reflections to consider or a single number of reflections.  For tropospheric paths at notable distance, it might be useful to start the eigenray search at a larger number of reflections and separate search runs are useful for identifying middle- and upper atmospheric returns that might have just a few reflections and tropospheric paths that might require more (e.g., stratospheric returns with 1 or 2 reflections and tropospheric paths with 10 - 15 reflections).

  *Search parameters*: the initial search for eigenrays uses an adustable inclination step that can be controlled using :code:`--incl-step-min` and :code:`--incl-step-max` (defaults 0.001 and 0.1 degrees).  Once an estimated eigenray is identified, the azimuth deviation due to cross winds is checked against the :code:`--az-dev-lim` value to ensure the LM search will be well posed and stable  (defaults to 2 degrees).

  *LM parameters*: the Levenberg-Marquardt algorithm uses a damping coefficient to stabilize the search which can be modified using :code:`--damping` (:math:`\lambda` in the discussion of :ref:`physics`).  A tolerance for accepting an eigenray can be accessed through :code:`--tolerance`, which has a default value of 0.1 km (smaller than the typical infrasonic array aperture).

Eigenray results can be visualized using :code:`infraga plot eigenray --atmo-file ToyAtmo.met` and show the along-azimuth effective sound speed profile for reference, all identified eigenray paths, and some characteristic (default is launch inclination vs. travel time).  As seen below, the color of the ray path corresponds to the arrival characteristics so the individual paths can be compared.  In this case, a direct thermospheric path (blue) is identified with later arrival time than the pair of single-reflection stratospheric arrivals (orange, green).

  .. image:: _static/_images/eigenray1.png
      :width: 1200px
      :align: center


The :code:`--y-axis-option` parameter included in the azimuthal visualizations is available in this method as well.  Running the above visualization and adding :code:`--y-axis-option amplitude` produces the below result so one can see how the largest amplitude amplitude is the shallower angle stratospheric arrival (often termed the 'slow stratospheric phase' per Waxler et al. (2015)).  The earlier, steeper inclination angle stratospheric phase (often termed the 'fast stratospheric phase') is markedly lower amplitude, and the thermospheric phase arriving later is even more attenuated due to combination of geometric and thermo-viscous losses.

  .. image:: _static/_images/eigenray2.png
      :width: 1200px
      :align: center

**Waveform Calculations**

Waveform evolution along individual ray paths can be computed in the weakly non-linear (also termed the weak shock limit) using a Burgers equation method (summarized in the discussion of :ref:`physics` and in Lonzaga et al. (2015) and Blom & Waxler (2021)) through the :code:`wnl-wvfrm` option.  Such analysis requires a known waveform at some location, :math:`s_1`, along the ray path and computes the evolved waveform at some later location, :math:`s_2` further along the path.  Because such analysis focuses on the waveform evolution along a single ray path, a single inclination and azimuth angle combination are needed to identify the ray path of interest.

Several waveform options are included in infraGA/GeoAc including a blastwave/impulse (originally from Freidlander (1946) and generalized by Waxler et al. (2018), see :ref:`physics` discussion for details) and an "N-wave" typical of sonic booms.  Alternatively, one can specify a file containing an overpressure time series recorded at some location and use that to initialize the waveform evolution methods.  Lastly, though not included in the C/C++ interface, the Python CLI includes an option to specify an explosive yield in equivalent kg of TNT for which the Kinney & Graham (1985) scaling laws are used to initialize a blastwave.  

Using the slow stratospheric phase in the above eigenray examples (theta, phi = 4.1392114, -84.970504 degrees; 1 ground reflection) and considering a 10 ton (10e3 kg) eq. TNT source,

  .. code:: none

    infraga sph wnl-wvfrm --atmo-file ToyAtmo.met --src-lat 30.0 --src-lon -100.0 --bounces 1 --inclination 4.1392114 --azimuth -84.970504 --wvfrm-yield 10e3

  .. code:: none

      ###########################################
      ####        Running infraga-sph        ####
      ####     Weakly Non-Linear Waveform    ####
      ###########################################

    Interpolating atmosphere data in 'ToyAtmo.met' using format 'zTuvdp'...
      Propagation region limits:
        latitude = -90, 90
        longitude = -180, 180
        altitude = 0, 139.9


    Parameter summary:
      inclination: 4.13921
      azimuth: -84.9705
      bounces: 1
      source location (lat, lon, alt): 30, -100, 0
      ground elevation: 0
      frequency: 0.1
      S&B atten coeff: 1
      waveform option: impulse
      waveform p0: 2432.82
      waveform t0: 0.0890579
      waveform alpha: 0.01
      waveform reference location: 0.754052

    Calculating ray path geometry and weakly non-linear waveform evolution...
      Defining waveform from built in impulse with shaping parameter = 0.01.
      Generating waveform at ray length: 0	altitude: 0.0579589	scaled non-linearity factor: 0.0458741
      Generating waveform at ray length: 10	altitude: 1.32843	scaled non-linearity factor: 0.00360202
      ...
      Generating waveform at ray length: 110	altitude: 37.596	scaled non-linearity factor: 0.026272
      Caustic encountered.
      Generating waveform at ray length: 120	altitude: 36.9643	scaled non-linearity factor: 0.0178894
      ...
      Generating waveform at ray length: 440	altitude: 0.347231	scaled non-linearity factor: 0.000299998

      Arrival summary:
        latitude [deg] = 30.249997
        longitude [deg] = -104.25002
        time [s] = 1412.4271
        celerity [km/s] = 0.29000158
        turning height [km] = 37.596484
        arrival inclination [deg] = 4.1765906
        back azimuth [deg] = 92.895927
        trans. coeff. [dB] = -40.225568
        absorption [dB] = -0.14986747


The peak overpressure and positive phase duration are defined from Kinney & Graham (1985) scaling laws computed internally at a reference distance defined by 35 :math:`\frac{\text{m}}{\text{kg}^3}` (this is the distance at which the shock velocity slows to within 1% of the sound speed and propagation becomes essentially linear/elastic).  Also, the ambient pressure and temperature are estimated from the sound speed and density so that explosive sources aloft can be specified.  The waveform evolution is computed along the ray path and any caustic encounters are noted and the expected phase shift applied.  The analysis returns an initial waveform file ('ToyAtmo.wvfrm_init.dat') that includes the near-source blastwave and the arrival waveform file ('ToyAtmo.wvfrm_out.dat').  

Currently there are no built-in visualization methods for individual waveform predictions.

**Combined Eigenray + Waveform Methods**

The C/C++ methods in infraGA/GeoAc for computing eigenrays and the weakly non-linear waveform methods are separated, but have been wrapped into a common Python method that identifies eigenrays connecting a specified source-receiver pair, checks for duplicate eigenray solutions (these can occur when ray paths are near critical launch angles that transition between waveguides), computes waveform predictions for each eigenray, and returns a merged set of results.  The parameter sets used for the eigenray and weakly non-linear waveform analyses are all available in this combined method along with a few additional parameters.  Waveform predictions for the above eigenray analysis can be completed by running,

  .. code-block:: none 

    infraga sph eig_wvfrm --atmo-file ToyAtmo.met --src-lat 30.0 --src-lon -100.0 --rcvr-lat 30.25 --rcvr-lon -104.25 --bnc-max 1 --keep-eig-results True --wvfrm-yield 10e3

The :code:`eig_wvfrm` methods run the eigenray and weakly non-linear methods and merge results into a pair of results files: 'ToyAtmo.eigenrays.dat' and 'ToyAtmo.wvfrms.dat'.  The eigenrays file includes the combined set of all identified ray paths.  The waveforms file has a variable number of columns containing the time (relative to the reference waveform) and a column containing the waveform contribution for each eigenray.  The header of the waveform results file includes much of the information from the simulation as well as the eigenray arrival information and the ordering of the arrival information matches that of the waveform contributions.

  .. code-block:: none

    # 'infraga sph eig_wvfrm' waveform results 
    #
    # 	profile: ToyAtmo.met
    # 	source location (lat, lon, alt): 30.0, -100.0, 0.0
    # 	receiver location (lat, lon, alt): 30.25, -104.25, 0.0
    # 	inclination range: 0.5, 45.0
    # 	inclination step max: 0.1
    # 	bounces: 0, 1
    # 	ground elevation: 0.0
    # 	damping: 0.001
    #   range max: 2500.0 
    #
    #    waveform reference distance: 0.7540521415111592
    #    waveform length: 0.05 None
    #    waveform source yield: 10e3
    #
    # Eigenray arrivals:
    # incl [deg]	az [deg]	n_b	lat_0 [deg]	lon_0 [deg]	time [s]	cel [km/s]	turning ht [km]	inclination [deg]	back azimuth [deg]	trans. coeff. [dB]	absorption [dB]
    # 35.831669 -84.874356 0.0 30.249983 -104.24998 1503.657 0.27240416 130.22471 35.84037 92.992416 -54.837386 -23.705732
    # 4.1392114 -84.970504 1.0 30.249997 -104.25002 1412.4271 0.29000158 37.596484 4.1765906 92.895927 -40.225568 -0.14986747
    # 31.703267 -84.583914 1.0 30.249995 -104.24994 1392.5105 0.29414393 48.989872 31.730785 93.283741 -59.050528 -0.29681437

    # t [s]	p1 [Pa] 	p2 [Pa] ...
    1334.9252	0.0	0.0	0.0	
    1334.9392	0.0	0.0	1.5782544e-12	
    1334.9532000000002	0.0	0.0	6.223907061041175e-12
    ...

The :code:`--keep-eig-results` option preserves the individual eigenray path files from the eigenray stage (useful if you want to run the eigenray visualization to see each arrival's characteristics) as well as the arrivals file itself.  It should be noted that keeping the arrivals file enables rapid calculation of new waveform predictions if the source is modified.  That is, if the above simulation has been run and one wants to see how the results change for a 100 ton eq. TNT source, one can simply modify the :code:`--wvfrm-yield` parameter and the Python methods will check whether eigenray results are already present and skip that portion of the analysis if results are found.

The combined eigenray and waveform predictions can be visualized using similar syntax to the azimuthal and eigenray visualizations :code:`infraga plot eig_wvfrms --atmo-file ToyAtmo.met`,

  .. image:: _static/_images/eig_wvfrm1.png
      :width: 600px
      :align: center

The resulting figures shows the eigenray paths, predicted overpressure waveform, and arrival characteristics using a matching travel time axis.  In this example, the stratospheric pair is easily identified in the early potion of the waveform and can be seen more clearly by zooming (again, click/drag the matplotlib window) as shown below. 

  .. image:: _static/_images/eig_wvfrm2.png
      :width: 600px
      :align: center

The later arrival from the upper atmosphere can also be more clearly visualized using the interactive matplotlib window as shown below.  In this case the relatively low amplitude of the source and the high turning height of the propagation path results in an overly small amplitude for this waveform contribution.

  .. image:: _static/_images/eig_wvfrm3.png
      :width: 600px
      :align: center

Lastly, similar to the azimuthal and eigenray visualization methods, the lower panel can be modified by specifying the y-axis and colormap variables.  Visualization of the trace velocity and amplitude information can be achieved as,

  .. code-block:: none
    
    infraga plot eig_wvfrms --atmo-file ToyAtmo.met --y-axis-option trace-velocity --cmap-option amplitude


  .. image:: _static/_images/eig_wvfrm4.png
      :width: 800px
      :align: center


**Supersonic Sources**

The set of ray paths emitted by a supersonic source as defined by the geometry of the Mach cone can be computed as discussed in Blom et al. (2024).  The Mach cone ray paths are computed for a single source location by specifying a source location, orientation (attack and azimuth angles), and Mach number (velocity relative to the ambient sound speed).  Usage of this method for a single point source can be run using the below syntax:

  .. code-block:: none

    infraga-sph -mach_cone G2S_example.met output_id=temp/t0_110.0 src_mach=2.63 src_attack=-13.12 src_az=90.45 cone_resol=1.0

In practice, computation of the infrasonic signals produced by a supersonic source's Mach cone requires computation of the above set of Mach cone rays at each point along a trajectory and then combination of all predicted paths with appropriate time delay.  The python interface for *infraGA/GeoAc* includes a method that will ingest atmospheric data as well as a trajectory file (containing time, latitude, longitude, and altitude).  The method computes the Mach number and source orientation information from the trajectory and steps through all points on the trajectory.  Running the prediction of rays for a supersonic source requires an atmospheric file (or range dependent grid) and a trajectory file:

  .. code-block:: none
    
    infraga sph supersonic --atmo-file G2S_example.met --trajectory trajectories/ballistic_traj.dat --traj-step 5 --cpu-cnt 4 --cleanup False --output-id ballistic

The :code:`--traj-step` allows one to skip through high resolution trajectory information (the step is that used in *numpy* indexing notation, :code:`traj_data[::k]`).  The :code:`--cleanup` flag can be used to keep the predicted ray path data for individual trajectory points.  In general, only arrival information is stored when running the method unless the :code:`--write-rays` option is turned on.  While ray paths are being computed, the trajectory information is displayed in a window with points showing where along the trajectory the current computations is located.

  .. image:: _static/_images/supersonic_trajectory.png
    :width: 800px
    :align: center 

Once the computation is complete, the arrival information can be visualized as with other ray tracing simulations:


  .. code-block:: none
    
    infraga plot map --arrivals ballistic.arrivals.dat


  .. image:: _static/_images/ballistic_arrivals.png
    :width: 600px
    :align: center 

More information about simulating infrasound from supersonic sources can be found in Blom et al. (2024) and citations therein.
