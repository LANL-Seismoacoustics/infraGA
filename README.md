# infraGA/GeoAc
InfraGA/GeoAc is a tool for modeling propagation of infrasound in the limit of geometric acoustics.  It includes methods for 2D and 3D Cartesian as well as spherical coordinate systems.  Eigenray methods are included to identify propagation paths between specific source-receiver geometries in any of the available coordiante systems.  Options are included to allow for range dependence in the atmosphere, realistic terrain, and accelerate methods via parallelization.

## Authorship
InfraGA/GeoAc is authored by Dr. Philip Blom

## Build
The standard methods can be compiled using a simple "make" command.
Accelerated methods are available using OpenMPI and can be compiled using "make accel".
