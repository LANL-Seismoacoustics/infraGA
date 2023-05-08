
=============================
Overview
=============================


infraGA/GeoAc is a set of tools for modeling the propagation of **infra**\ sound signals through the atmosphere in the limit of **Geo**\ metric **Ac**\ oustics.  The propagation methods are implemented in C/C++ with Python wrappers, utilities, and visualization tools.  infraGA/GeoAc includes tools to model the propagation of infrasonic signals through Cartesian geometry in 2D using the effective sound speed approximation or in 3D using a Cartesian moving medium atmosphere as well as using a spherical atmospheric layer surrounding the globe.  Algorithms include general point-source propagation, identification of eigenray paths connecting specific source-receiver geometries, as well as weakly non-linear waveform computations along ray paths.

InfraGA/GeoAc utilizes a set of auxiliary parameters to solve the Transport equation as discussed in Blom & Waxler (2012) that enables calculation of geometric spreading losses and leading order amplitude estimation.  These auxiliary parameters are further leveraged in the eigenray analysis as they provide an efficient means of computing the Jacobian for the ray path arrival location (see Blom & Waxler, 2017).  Weakly non-linear waveform evolution is computed using a Burgers' equation method as detailed in Lonzaga et al. (2015) and Blom & Waxler (2021).  Recent R&D includes the use of non-flat ground in modeling propagation (see Blom 2020) and modeling signals produced by the Mach cone of supersonic sources (in development).


**License**

Copyright (c) 2014, Triad National Security, LLC

All rights reserved.

Copyright 2014. Triad National Security, LLC. This software was produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, reproduce, and distribute this software. NEITHER THE GOVERNMENT NOR TRIAD NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE. If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL.

Additionally, permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

