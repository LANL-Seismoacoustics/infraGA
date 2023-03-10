.. _quickstart:

=====================================
Quickstart
=====================================

****************************
Command Line Interface (CLI) 
****************************

Most of InfraPy's analysis methods are accessible through a command line interface (CLI) with parameters specified either via command line flags or a configuration file.  Waveform data can be ingested from local files (eg., SAC or similar format that can be ingested via :code:`obspy.core.read`) or downloaded from FDSN clients via :code:`obspy.clients.fdsn`.  Array- and network-level analyses can be performed from the command line enabling a full pipeline of analysis from beamforming/detection to event identification and localization.  Visualization methods are also included to quickly interrogate analysis results.  The Quickstart summarized here steps through these various CLI methods and demonstrates the usage of InfraPy from the command line.

--------------------
Array-Level Analyses
--------------------

- The beamforming methods in InfraPy can be run via the :code:`run_fk` CLI option.  For a local data source such as the included SAC files in the data directory, this is simply,

    .. code-block:: bash

        infrapy run_fk --local-wvfrms 'data/YJ.BRP*.SAC'

    Note that the data path must be in quotes in order to properly parsed and that this Quickstart assumes you are in the infrapy/examples directory (if you are getting an error that the waveform data isn't found, make sure you're in the correct directory).  As the methods are run, data and algorithm parameters are summarized and a progress bar shows how much of the data has been analyzed:

    .. code-block:: none

        #####################################
        ##                                 ##
        ##             InfraPy             ##
        ##    Beamforming (fk) Analysis    ##
        ##                                 ##
        #####################################

        Data parameters:
          local_wvfrms: data/YJ.BRP*.SAC
          local_latlon: None
          local_fk_label: None

        Algorithm parameters:
          freq_min: 0.5
          freq_max: 5.0
          back_az_min: -180.0
          back_az_max: 180.0
          back_az_step: 2.0
          trace_vel_min: 300.0
          trace_vel_max: 600.0
          trace_vel_step: 2.5
          method: bartlett
          signal_start: None
          signal_end: None
          window_len: 10.0
          sub_window_len: None
          window_step: 5.0

        Loading local data from data/YJ.BRP*.SAC

        Data summary:
        YJ.BRP1..EDF	2012-04-09T18:00:00.008300Z - 2012-04-09T18:19:59.998300Z
        YJ.BRP2..EDF	2012-04-09T18:00:00.008300Z - 2012-04-09T18:19:59.998300Z
        YJ.BRP3..EDF	2012-04-09T18:00:00.008300Z - 2012-04-09T18:19:59.998300Z
        YJ.BRP4..EDF	2012-04-09T18:00:00.008300Z - 2012-04-09T18:19:59.998300Z

        Running fk analysis...
	        Progress: [>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]

        Writing results into data/YJ.BRP_2012.04.09_18.00.00-18.19.59.fk_results.dat
