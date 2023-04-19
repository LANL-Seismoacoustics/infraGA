.. _parameters:

=====================================
Parameters
=====================================

******************
General Parameters
******************

**Common Parameters (all geometries, all options)**

All infraga/GeoAc simulations have a common set of parameters that can be modified as summarized in the below table.  Note that there is a syntax difference when specifying a parameter value in the C/C++ command line interface (e.g., :code:`freq=0.1`) compared with the Python interface that uses the Click library (eg., :code:`--freq 0.1`).  Somewhat confusingly, the way in which the Python config file parsing works requires parameters defined in a configuration file to use the C/C++ syntax (underscores replacing hyphens).  Most of these parameters are relatively self explanatory or detailed in the :ref:`quickstart` or :ref:`advanced` sections.  All parameters in this section are defined in the :code:`[GENERAL]` section of configuration file.

+--------------------+----------------------+------------+----------------------+
| **Parameter**                             | **Info**                          |
+--------------------+----------------------+------------+----------------------+
| **C/C++ CLI**      | **Python CLI**       | **Units**  | **Default Value**    | 
+--------------------+----------------------+------------+----------------------+
| :code:`freq`       | :code:`--freq`       | Hz         | 0.1                  |
+--------------------+----------------------+------------+----------------------+
| :code:`abs_coeff`  | :code:`--abs-coeff`  | Scalar     | 1.0                  |
+--------------------+----------------------+------------+----------------------+
| :code:`z_grnd`     | :code:`--z-grnd`     | km         | 0.0                  |
+--------------------+----------------------+------------+----------------------+
| :code:`write_atmo` | :code:`--write-atmo` | True/False | False                |
+--------------------+----------------------+------------+----------------------+
| :code:`prof_format`| :code:`--prof-format`| string     | 'zTuvdp'             |
+--------------------+----------------------+------------+----------------------+
| :code:`output_id`  | :code:`--output-id`  | string     | From atmosphere file |
+--------------------+----------------------+------------+----------------------+
| :code:`calc_amp`   | :code:`--calc-amp`   | True/False | True                 |
+--------------------+----------------------+------------+----------------------+
| :code:`max_alt`    | :code:`--max-alt`    | km         | From atmosphere file |
+--------------------+----------------------+------------+----------------------+
| :code:`max_rng`    | :code:`--max-rng`    | km         | 1000                 |
+--------------------+----------------------+------------+----------------------+
| :code:`min_ds`     | :code:`--min-ds`     | km         | 0.001                |
+--------------------+----------------------+------------+----------------------+
| :code:`max_ds`     | :code:`--max-ds`     | km         | 0.05                 |
+--------------------+----------------------+------------+----------------------+
| :code:`max_s`      | :code:`--max-s`      | km         | 1000                 |
+--------------------+----------------------+------------+----------------------+
| :code:`topo_file`  | :code:`--topo-file`  | string     | See :ref:`advanced`  |
+--------------------+----------------------+------------+----------------------+
| **-**              | :code:`--config-file`| string     | None                 |
+--------------------+----------------------+------------+----------------------+

**3D Cartesian Region Boundaries**

For all 3D Cartessian anlayses, the :math:`x` (east/west) and :math:`y` (north/south) edges of the propagation domain can be defined:

+--------------------+----------------------+------------+----------------------+
| **Parameter**                             | **Info**                          |
+--------------------+----------------------+------------+----------------------+
| **C/C++ CLI**      | **Python CLI**       | **Units**  | **Default Value**    | 
+--------------------+----------------------+------------+----------------------+
| :code:`min_x`      | :code:`--min-x`      | km         | -1000                |
+--------------------+----------------------+------------+----------------------+
| :code:`max_x`      | :code:`--max-x`      | km         | 1000                 |
+--------------------+----------------------+------------+----------------------+
| :code:`min_y`      | :code:`--min-y`      | km         | -1000                |
+--------------------+----------------------+------------+----------------------+
| :code:`max_y`      | :code:`--max-y`      | km         | 1000                 |
+--------------------+----------------------+------------+----------------------+

**Spherical Geometry Region Boundaries**

In the current implementation of the spherical atmospheric layer propagation geometry, wrapping of ray paths across the global poles is not enabled so that ray paths are terminated if they extend beyond latitudes of :math:`\pm90^\circ`.  Wrapping of ray paths at the :math:`\pm180^\circ` longitude boundary is included in analysis so that keeping the below bounds allows ray paths to wrap through this boundary.  It should be noted that some of the visualization methods included don't work well for data that crosses this boundary.  For a smaller region of the globe, the latitude and longitude bounds can be defined as below.

+--------------------+----------------------+------------+----------------------+
| **Parameter**                             | **Info**                          |
+--------------------+----------------------+------------+----------------------+
| **C/C++ CLI**      | **Python CLI**       | **Units**  | **Default Value**    | 
+--------------------+----------------------+------------+----------------------+
| :code:`min_lat`    | :code:`--min-lat`    | degress    | -90                  |
+--------------------+----------------------+------------+----------------------+
| :code:`max_lat`    | :code:`--max-lat`    | degress    | 90                   |
+--------------------+----------------------+------------+----------------------+
| :code:`min_lon`    | :code:`--min-lon`    | degress    | -180                 |
+--------------------+----------------------+------------+----------------------+
| :code:`max_lon`    | :code:`--max-lon`    | degress    | 180                  |
+--------------------+----------------------+------------+----------------------+

**Weakly Non-Linear Waveform Evolution**

Weakly non-linear waveform methods require a single ray path and initialization of the waveform at some location near the source in order to compute the evolution of the waveform along the ray path.  The initial waveform can be specified from a file or through the built in options ('impulse', 'u-wave', or 'n-wave').  The Python interface includes an implementation of the Kinney & Graham scaling laws that enables use of an explosive yield (kg eq. TNT) to initialize the waveform.  The script checks the source altitude and atmosphere file to define the ambient temperature and pressure information and computes the peak overpressure and positive phase duration for the impulse (blastwave) near the source.  See :ref:`advanced` for a full discussion of waveform calculation usage.  The parameters below are specified in the :code:`[WAVEFORM]` section of configuration files.

+-------------------------------------------------+------------------------------------+
| **Parameter**                                   | **Info**                           |
+-----------------------+-------------------------+-------------+----------------------+
| **C/C++ CLI**         | **Python CLI**          | **Units**   | **Default Value**    | 
+-----------------------+-------------------------+-------------+----------------------+
| :code:`inclination`   | :code:`--inclination`   | degress     | 15.0                 |
+-----------------------+-------------------------+-------------+----------------------+
| :code:`azimuth`       | :code:`--azimuth`       | degress     | -90.0                |
+-----------------------+-------------------------+-------------+----------------------+
| :code:`bounces`       | :code:`--bounces`       | integer     | 10                   |
+-----------------------+-------------------------+-------------+----------------------+
| :code:`write_ray`     | :code:`--write-ray`     | True/False  | False                |
+-----------------------+-------------------------+-------------+----------------------+
| :code:`wvfrm_file`    | :code:`--wvfrm-file`    | string      | None                 |
+-----------------------+-------------------------+-------------+----------------------+
| :code:`wvfrm_opt`     | :code:`--wvfrm-opt`     | string      | 'impulse'            |
+-----------------------+-------------------------+-------------+----------------------+
| :code:`wvfrm_p0`      | :code:`--wvfrm-p0`      | Pa          | 10.0                 |
+-----------------------+-------------------------+-------------+----------------------+
| :code:`wvfrm_t0`      | :code:`--wvfrm-p0`      | seconds     | 1.0                  |
+-----------------------+-------------------------+-------------+----------------------+
| :code:`wvfrm_alpha`   | :code:`--wvfrm-p0`      | scalar      | 1.0                  |
+-----------------------+-------------------------+-------------+----------------------+
| :code:`wvfrm_ref`     | :code:`--wvfrm-p0`      | km          | 1.0                  |
+-----------------------+-------------------------+-------------+----------------------+
| :code:`wvfrm_out_step`| :code:`--wvfrm-out-step`| km          | None                 |
+-----------------------+-------------------------+-------------+----------------------+
| :code:`wvfrm_ds`      | :code:`--wvfrm-ds`      | km          | 0.5                  |
+-----------------------+-------------------------+-------------+----------------------+
| :code:`wvfrm_len`     | :code:`--wvfrm-len`     | int         | :math:`2^{13}`       |
+-----------------------+-------------------------+-------------+----------------------+
| **-**                 | :code:`--wvfrm-yield`   | kg (eq. TNT)| None                 |
+-----------------------+-------------------------+-------------+----------------------+

*************************************
Two-Dimensional Simulation Parameters
*************************************

Simulation of ray paths in an azimuthal plane (range vs. altitude) using the effective sound speed approximation is parameterized by a set of inclination angles, a single azimuth, a maximum number of possible ground bounces (reflections) and the source altitude (relative to sea level).  Waveform analysis in 2D requires only the source altitude in addition to the above common parameters.  In the below tables, the Point Source Propagation parameters are defined in the :code:`[PROP]` section of configuration files and the waveform parameters are defined using the :code:`[WAVEFORM]` header.


**Point Source Propagation**

+--------------------+----------------------+------------+----------------------+
| **Parameter**                             | **Info**                          |
+--------------------+----------------------+------------+----------------------+
| **C/C++ CLI**      | **Python CLI**       | **Units**  | **Default Value**    | 
+--------------------+----------------------+------------+----------------------+
| :code:`incl_min`   | :code:`--incl-min`   | degress    | 0.5                  |
+--------------------+----------------------+------------+----------------------+
| :code:`incl_max`   | :code:`--incl-max`   | degress    | 45.0                 |
+--------------------+----------------------+------------+----------------------+
| :code:`incl_step`  | :code:`--incl-step`  | degress    | 0.5                  |
+--------------------+----------------------+------------+----------------------+
| :code:`inclination`| :code:`--inclination`| degress    | None                 |
+--------------------+----------------------+------------+----------------------+
| :code:`azimuth`    | :code:`--azimuth`    | degress    | -90.0                |
+--------------------+----------------------+------------+----------------------+
| :code:`bounces`    | :code:`--bounces`    | integer    | 10                   |
+--------------------+----------------------+------------+----------------------+
| :code:`src_alt`    | :code:`--src-alt`    | km         | 0.0                  |
+--------------------+----------------------+------------+----------------------+

**Weakly Non-Linear Waveform Evolution**

+-------------------------------------------------+-----------------------------------+
| **Parameter**                                   | **Info**                          |
+-----------------------+-------------------------+------------+----------------------+
| **C/C++ CLI**         | **Python CLI**          | **Units**  | **Default Value**    | 
+-----------------------+-------------------------+------------+----------------------+
| :code:`src_alt`       | :code:`--src-alt`       | km         | 0.0                  |
+-----------------------+-------------------------+------------+----------------------+

***************************************
Three-Dimensional Simulation Parameters
***************************************

Computation of ray paths in 3D (Cartesian) geometry requires both a set of inclination angles as well as azimuth angles.  The 'azimuth' and 'inclination' parameters define single values of the respective angles for simulations focused on a single azimuthal direction or considering propagation at all compass directions from a single radiating inclination, respectively.  The source location for 3D Cartesian simulations requires an :math:`x` (east/west), :math:`y` (north/south), and altitude.  Because multi-azimuth simulations produce relatively large volumes of ray path data, an option is included to prevent that information from being written to file and instead only return the arrival information where ray paths intercept the ground surface.  In the case of eigenray analysis, both the source and receiver locations are defined in terms of their :math:`x` (east/west), :math:`y` (north/south) locations.  The receiver location is assumed to be on the ground surface, but the source location can be on the surface or aloft.  The various other eigenray parmaeters are summarized in :ref:`advanced`.  In the below tables, the Point Source Propagation parameters are defined in the :code:`[PROP]` section of configuration files, Eigenray Analysis parameters using the :code:`[EiGENRAY]` header, and the waveform parameters are defined using the :code:`[WAVEFORM]` header.

**Point Source Propagation**

+--------------------+----------------------+------------+----------------------+
| **Parameter**                             | **Info**                          |
+--------------------+----------------------+------------+----------------------+
| **C/C++ CLI**      | **Python CLI**       | **Units**  | **Default Value**    | 
+--------------------+----------------------+------------+----------------------+
| :code:`incl_min`   | :code:`--incl-min`   | degress    | 0.5                  |
+--------------------+----------------------+------------+----------------------+
| :code:`incl_max`   | :code:`--incl-max`   | degress    | 45.0                 |
+--------------------+----------------------+------------+----------------------+
| :code:`incl_step`  | :code:`--incl-step`  | degress    | 0.5                  |
+--------------------+----------------------+------------+----------------------+
| :code:`inclination`| :code:`--inclination`| degress    | None                 |
+--------------------+----------------------+------------+----------------------+
| :code:`az_min`     | :code:`--az-min`     | degress    | -90.0                |
+--------------------+----------------------+------------+----------------------+
| :code:`az_max`     | :code:`--az-max`     | degress    | -90.0                |
+--------------------+----------------------+------------+----------------------+
| :code:`az_step`    | :code:`--az-step`    | degress    | 1.0                  |
+--------------------+----------------------+------------+----------------------+
| :code:`azimuth`    | :code:`--azimuth`    | degress    | -90.0                |
+--------------------+----------------------+------------+----------------------+
| :code:`bounces`    | :code:`--bounces`    | integer    | 10                   |
+--------------------+----------------------+------------+----------------------+
| :code:`src_x`      | :code:`--src-x`      | km         | 0.0                  |
+--------------------+----------------------+------------+----------------------+
| :code:`src_y`      | :code:`--src-y`      | km         | 0.0                  |
+--------------------+----------------------+------------+----------------------+
| :code:`src_alt`    | :code:`--src-alt`    | km         | 0.0                  |
+--------------------+----------------------+------------+----------------------+
| :code:`turn_ht_min`| :code:`--turn-ht-min`| km         | 0.2                  |
+--------------------+----------------------+------------+----------------------+
| :code:`write_rays` | :code:`--write-rays` | True/False | True                 |
+--------------------+----------------------+------------+----------------------+
| :code:`write_topo` | :code:`--write-topo` | True/False | False                |
+--------------------+----------------------+------------+----------------------+

**Eigenray Analysis**

+-----------------------------------------------+-----------------------------------+
| **Parameter**                                 | **Info**                          |
+----------------------+------------------------+------------+----------------------+
| **C/C++ CLI**        | **Python CLI**         | **Units**  | **Default Value**    | 
+----------------------+------------------------+------------+----------------------+
| :code:`incl_min`     | :code:`--incl-min`     | degress    | 0.5                  |
+----------------------+------------------------+------------+----------------------+
| :code:`incl_max`     | :code:`--incl-max`     | degress    | 45.0                 |
+----------------------+------------------------+------------+----------------------+
| :code:`bnc_min`      | :code:`--bnc-min`      | integer    | 0                    |
+----------------------+------------------------+------------+----------------------+
| :code:`bnc_max`      | :code:`--bnc-max`      | integer    | 0                    |
+----------------------+------------------------+------------+----------------------+
| :code:`bounces`      | :code:`--bounces`      | integer    | 0                    |
+----------------------+------------------------+------------+----------------------+
| :code:`src_x`        | :code:`--src-x`        | km         | 0.0                  |
+----------------------+------------------------+------------+----------------------+
| :code:`src_y`        | :code:`--src-y`        | km         | 0.0                  |
+----------------------+------------------------+------------+----------------------+
| :code:`src_alt`      | :code:`--src-alt`      | km         | 0.0                  |
+----------------------+------------------------+------------+----------------------+
| :code:`rcvr_x`       | :code:`--rcvr-x`       | km         | 250.0                |
+----------------------+------------------------+------------+----------------------+
| :code:`rcvr_y`       | :code:`--rcvr-y`       | km         | 0.0                  |
+----------------------+------------------------+------------+----------------------+
| :code:`verbose`      | :code:`--verbose`      | True/False | False                |
+----------------------+------------------------+------------+----------------------+
| :code:`iterations`   | :code:`--iterations`   | integer    | 25                   |
+----------------------+------------------------+------------+----------------------+
| :code:`damping`      | :code:`--damping`      | scalar     | 1.0e-3               |
+----------------------+------------------------+------------+----------------------+
| :code:`tolerance`    | :code:`--tolerance`    | km         | 0.1                  |
+----------------------+------------------------+------------+----------------------+
| :code:`az_dev_lim`   | :code:`--az-dev-lim`   | degrees    | 2.0                  |
+----------------------+------------------------+------------+----------------------+
| :code:`incl_step_min`| :code:`--incl-step-min`| degrees    | 0.001                |
+----------------------+------------------------+------------+----------------------+
| :code:`incl_step_max`| :code:`--incl-step-max`| degrees    | 0.1                  |
+----------------------+------------------------+------------+----------------------+

**Weakly Non-Linear Waveform Evolution**

+--------------------+----------------------+------------+----------------------+
| **Parameter**                             | **Info**                          |
+--------------------+----------------------+------------+----------------------+
| **C/C++ CLI**      | **Python CLI**       | **Units**  | **Default Value**    | 
+--------------------+----------------------+------------+----------------------+
| :code:`src_x`      | :code:`--src-x`      | km         | 0.0                  |
+--------------------+----------------------+------------+----------------------+
| :code:`src_y`      | :code:`--src-y`      | km         | 0.0                  |
+--------------------+----------------------+------------+----------------------+
| :code:`src_alt`    | :code:`--src-alt`    | km         | 0.0                  |
+--------------------+----------------------+------------+----------------------+


************************************************
Spherical Atmosphere Layer Simulation Parameters
************************************************

Ray path calculation in a spherical atmospheric layer geometry again utilizes sets of inclination and azimuthal angles and single values can be specified using 'inclination' or 'azimuth', respectively.  The location of the source (and receiver for eigenray analysis) is defined by the latitude and longitude on the globe.  In most other respects, the parameter set for the spherical geometry methods is identical to that of the 3D Cartesian.  Similar to the 3D Cartesian configuration file header, the Point Source Propagation parameters are defined in the :code:`[PROP]` section of configuration files, Eigenray Analysis parameters using the :code:`[EiGENRAY]` header, and the waveform parameters are defined using the :code:`[WAVEFORM]` header.

**Point Source Propagation**

+--------------------+----------------------+------------+----------------------+
| **Parameter**                             | **Info**                          |
+--------------------+----------------------+------------+----------------------+
| **C/C++ CLI**      | **Python CLI**       | **Units**  | **Default Value**    | 
+--------------------+----------------------+------------+----------------------+
| :code:`incl_min`   | :code:`--incl-min`   | degress    | 0.5                  |
+--------------------+----------------------+------------+----------------------+
| :code:`incl_max`   | :code:`--incl-max`   | degress    | 45.0                 |
+--------------------+----------------------+------------+----------------------+
| :code:`incl_step`  | :code:`--incl-step`  | degress    | 0.5                  |
+--------------------+----------------------+------------+----------------------+
| :code:`inclination`| :code:`--inclination`| degress    | None                 |
+--------------------+----------------------+------------+----------------------+
| :code:`az_min`     | :code:`--az-min`     | degress    | -90.0                |
+--------------------+----------------------+------------+----------------------+
| :code:`az_max`     | :code:`--az-max`     | degress    | -90.0                |
+--------------------+----------------------+------------+----------------------+
| :code:`az_step`    | :code:`--az-step`    | degress    | 1.0                  |
+--------------------+----------------------+------------+----------------------+
| :code:`azimuth`    | :code:`--azimuth`    | degress    | -90.0                |
+--------------------+----------------------+------------+----------------------+
| :code:`bounces`    | :code:`--bounces`    | integer    | 10                   |
+--------------------+----------------------+------------+----------------------+
| :code:`src_lat`    | :code:`--src-lat`    | km         | 30.0                 |
+--------------------+----------------------+------------+----------------------+
| :code:`src_lon`    | :code:`--src-lon`    | km         | 0.0                  |
+--------------------+----------------------+------------+----------------------+
| :code:`src_alt`    | :code:`--src-alt`    | km         | 0.0                  |
+--------------------+----------------------+------------+----------------------+
| :code:`turn_ht_min`| :code:`--turn-ht-min`| km         | 0.2                  |
+--------------------+----------------------+------------+----------------------+
| :code:`write_rays` | :code:`--write-rays` | True/False | True                 |
+--------------------+----------------------+------------+----------------------+
| :code:`write_topo` | :code:`--write-topo` | True/False | False                |
+--------------------+----------------------+------------+----------------------+

**Eigenray Analysis**

+-----------------------------------------------+-----------------------------------+
| **Parameter**                                 | **Info**                          |
+----------------------+------------------------+------------+----------------------+
| **C/C++ CLI**        | **Python CLI**         | **Units**  | **Default Value**    | 
+----------------------+------------------------+------------+----------------------+
| :code:`incl_min`     | :code:`--incl-min`     | degress    | 0.5                  |
+----------------------+------------------------+------------+----------------------+
| :code:`incl_max`     | :code:`--incl-max`     | degress    | 45.0                 |
+----------------------+------------------------+------------+----------------------+
| :code:`bnc_min`      | :code:`--bnc-min`      | integer    | 0                    |
+----------------------+------------------------+------------+----------------------+
| :code:`bnc_max`      | :code:`--bnc-max`      | integer    | 0                    |
+----------------------+------------------------+------------+----------------------+
| :code:`bounces`      | :code:`--bounces`      | integer    | 0                    |
+----------------------+------------------------+------------+----------------------+
| :code:`src_lat`      | :code:`--src-lat`      | degrees    | 30.0                 |
+----------------------+------------------------+------------+----------------------+
| :code:`src_lon`      | :code:`--src-lon`      | degrees    | 0.0                  |
+----------------------+------------------------+------------+----------------------+
| :code:`src_alt`      | :code:`--src-alt`      | km         | 0.0                  |
+----------------------+------------------------+------------+----------------------+
| :code:`rcvr_lat`     | :code:`--rcvr-lat`     | degrees    | 30.0                 |
+----------------------+------------------------+------------+----------------------+
| :code:`rcvr_lon`     | :code:`--rcvr-lon`     | degrees    | 2.5                  |
+----------------------+------------------------+------------+----------------------+
| :code:`verbose`      | :code:`--verbose`      | True/False | False                |
+----------------------+------------------------+------------+----------------------+
| :code:`iterations`   | :code:`--iterations`   | integer    | 25                   |
+----------------------+------------------------+------------+----------------------+
| :code:`damping`      | :code:`--damping`      | scalar     | 1.0e-3               |
+----------------------+------------------------+------------+----------------------+
| :code:`tolerance`    | :code:`--tolerance`    | km         | 0.1                  |
+----------------------+------------------------+------------+----------------------+
| :code:`az_dev_lim`   | :code:`--az-dev-lim`   | degrees    | 2.0                  |
+----------------------+------------------------+------------+----------------------+
| :code:`incl_step_min`| :code:`--incl-step-min`| degrees    | 0.001                |
+----------------------+------------------------+------------+----------------------+
| :code:`incl_step_max`| :code:`--incl-step-max`| degrees    | 0.1                  |
+----------------------+------------------------+------------+----------------------+


**Weakly Non-Linear Waveform Evolution**


+--------------------+----------------------+------------+----------------------+
| **Parameter**                             | **Info**                          |
+--------------------+----------------------+------------+----------------------+
| **C/C++ CLI**      | **Python CLI**       | **Units**  | **Default Value**    | 
+--------------------+----------------------+------------+----------------------+
| :code:`src_lat`    | :code:`--src-lat`    | degrees    | 30.0                 |
+--------------------+----------------------+------------+----------------------+
| :code:`src_lon`    | :code:`--src-lon`    | degrees    | 0.0                  |
+--------------------+----------------------+------------+----------------------+
| :code:`src_alt`    | :code:`--src-alt`    | km         | 0.0                  |
+--------------------+----------------------+------------+----------------------+
