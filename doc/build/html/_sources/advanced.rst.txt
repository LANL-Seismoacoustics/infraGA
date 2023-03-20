.. _advanced:

=====================================
Advanced Usage
=====================================

Introduction...

****************************
Range Dependent Analysis 
****************************

Discussion of using a G2S specification grid...


  For even larger scale analyses (e.g., analysis of propagation paths from the Tonga volcanic eruption that circled the globe multiple times), one can define a maximum propagation range greater than that of the earth (~24,000 km) and use the :code:`--max-tm` parameter to define a maximum propagation time (defined in hours).

    .. code:: none

      infraga sph prop --atmo-file ToyAtmo.met --max-rng 25000 --max-tm 24.0 --bounces 1000

  It should be noted that for such large scales, the spatial variations in the atmosphere become significant and the range dependent methods detailed in the :ref:`advanced` discussion are needed.

  
****************************
Eigenray Analysis 
****************************

Discussion of eigenray methods...

****************************
Waveform Calculations 
****************************

Discussion of waveform methods...

************************************
Combined Eigenray + Waveform Methods
************************************

Discussion of eig_wvfrm methods...



****************************
Including Reaslistic Terrain 
****************************

Discussion of how to include terrain...

