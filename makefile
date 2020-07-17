CC= g++
mpiCC= mpic++
CFLAGS=  -Wno-write-strings -lfftw3

MAINS_DIR= src/
ATMO_DIR= src/atmo/
GEOAC_DIR= src/geoac/
UTILS_DIR= src/util/

INSTALL_DIR= /usr/local/bin

# Targets to build all of the standard methods and all of the accelerated methods

all: prep infraga-2d infraga-3d infraga-3d-rngdep infraga-sph infraga-sph-rngdep

accel: infraga-accel-3d infraga-accel-3d-rngdep infraga-accel-sph infraga-accel-sph-rngdep

# Declare which items don’t create anything

.PHONY: all accel clean clean_accel install install_accel uninstall uninstall_accel prep

# Make sure the bin directory exists

prep:
	mkdir -p bin

# Build the standard methods

infraga-2d: prep
	${CC}  ${ATMO_DIR}atmo_io.2d.cpp ${ATMO_DIR}atmo_state.2d.cpp ${UTILS_DIR}fileIO.cpp ${UTILS_DIR}interpolation.cpp ${ATMO_DIR}atmo_state.abs.cart.cpp ${GEOAC_DIR}geoac.params.2d.cpp ${GEOAC_DIR}geoac.eqset.2d.cpp ${UTILS_DIR}rk4solver.cpp ${GEOAC_DIR}geoac.interface.2d.cpp ${UTILS_DIR}waveforms.cpp ${MAINS_DIR}main.infraga.2d.cpp ${CFLAGS} -o bin/infraga-2d

infraga-3d: prep
	${CC}  ${ATMO_DIR}atmo_io.3d.strat.cpp ${ATMO_DIR}atmo_state.3d.strat.cpp ${UTILS_DIR}fileIO.cpp ${UTILS_DIR}interpolation.cpp ${ATMO_DIR}atmo_state.abs.cart.cpp ${GEOAC_DIR}geoac.params.3d.cpp ${GEOAC_DIR}geoac.eqset.3d.strat.cpp ${GEOAC_DIR}geoac.eigenray.3d.cpp ${UTILS_DIR}rk4solver.cpp ${GEOAC_DIR}geoac.interface.3d.cpp ${UTILS_DIR}waveforms.cpp ${MAINS_DIR}main.infraga.3d.cpp ${CFLAGS} -o bin/infraga-3d

infraga-3d-rngdep: prep
	${CC}  ${ATMO_DIR}atmo_io.3d.rngdep.cpp ${ATMO_DIR}atmo_state.3d.rngdep.cpp ${UTILS_DIR}fileIO.cpp ${UTILS_DIR}interpolation.cpp ${ATMO_DIR}atmo_state.abs.cart.cpp ${GEOAC_DIR}geoac.params.3d.cpp ${GEOAC_DIR}geoac.eqset.3d.rngdep.cpp ${GEOAC_DIR}geoac.eigenray.3d.cpp ${UTILS_DIR}rk4solver.cpp ${GEOAC_DIR}geoac.interface.3d.cpp ${UTILS_DIR}waveforms.cpp ${MAINS_DIR}main.infraga.3d.rngdep.cpp ${CFLAGS} -o bin/infraga-3d-rngdep

infraga-sph: prep
	${CC}  ${ATMO_DIR}atmo_io.sph.strat.cpp ${ATMO_DIR}atmo_state.sph.strat.cpp ${UTILS_DIR}globe.cpp ${UTILS_DIR}fileIO.cpp ${UTILS_DIR}interpolation.cpp ${ATMO_DIR}atmo_state.abs.sph.cpp ${GEOAC_DIR}geoac.params.sph.cpp ${GEOAC_DIR}geoac.eigenray.sph.cpp ${GEOAC_DIR}geoac.eqset.sph.strat.cpp ${UTILS_DIR}rk4solver.cpp ${GEOAC_DIR}geoac.interface.sph.cpp ${UTILS_DIR}waveforms.cpp ${MAINS_DIR}main.infraga.sph.cpp ${CFLAGS} -o bin/infraga-sph

infraga-sph-rngdep: prep
	${CC}  ${ATMO_DIR}atmo_io.sph.rngdep.cpp ${ATMO_DIR}atmo_state.sph.rngdep.cpp ${UTILS_DIR}globe.cpp ${UTILS_DIR}fileIO.cpp ${UTILS_DIR}interpolation.cpp ${ATMO_DIR}atmo_state.abs.sph.cpp ${GEOAC_DIR}geoac.params.sph.cpp ${GEOAC_DIR}geoac.eigenray.sph.cpp ${GEOAC_DIR}geoac.eqset.sph.rngdep.cpp ${UTILS_DIR}rk4solver.cpp ${GEOAC_DIR}geoac.interface.sph.cpp ${UTILS_DIR}waveforms.cpp ${MAINS_DIR}main.infraga.sph.rngdep.cpp ${CFLAGS} -o bin/infraga-sph-rngdep

# Build the accelerated methods (might require modifying epic above

infraga-accel-3d:
	${mpiCC}  ${ATMO_DIR}atmo_io.3d.strat.mpi.cpp ${ATMO_DIR}atmo_state.3d.strat.cpp ${UTILS_DIR}fileIO.cpp ${UTILS_DIR}interpolation.cpp ${ATMO_DIR}atmo_state.abs.cart.cpp ${GEOAC_DIR}geoac.params.3d.cpp ${GEOAC_DIR}geoac.eqset.3d.strat.cpp ${GEOAC_DIR}geoac.eigenray.3d.mpi.cpp ${UTILS_DIR}rk4solver.cpp ${GEOAC_DIR}geoac.interface.3d.cpp ${UTILS_DIR}waveforms.cpp ${UTILS_DIR}multitasking.cpp ${MAINS_DIR}main.infraga.3d.mpi.cpp ${CFLAGS} -o bin/infraga-accel-3d

infraga-accel-3d-rngdep:
	${mpiCC}  ${ATMO_DIR}atmo_io.3d.rngdep.mpi.cpp ${ATMO_DIR}atmo_state.3d.rngdep.cpp ${UTILS_DIR}fileIO.cpp ${UTILS_DIR}interpolation.cpp ${ATMO_DIR}atmo_state.abs.cart.cpp ${GEOAC_DIR}geoac.params.3d.cpp ${GEOAC_DIR}geoac.eqset.3d.rngdep.cpp ${GEOAC_DIR}geoac.eigenray.3d.mpi.cpp ${UTILS_DIR}rk4solver.cpp ${GEOAC_DIR}geoac.interface.3d.cpp ${UTILS_DIR}waveforms.cpp ${UTILS_DIR}multitasking.cpp ${MAINS_DIR}main.infraga.3d.rngdep.mpi.cpp ${CFLAGS} -o bin/infraga-accel-3d-rngdep

infraga-accel-sph:
	${mpiCC}  ${ATMO_DIR}atmo_io.sph.strat.mpi.cpp ${ATMO_DIR}atmo_state.sph.strat.cpp ${UTILS_DIR}globe.cpp ${UTILS_DIR}fileIO.cpp ${UTILS_DIR}interpolation.cpp ${ATMO_DIR}atmo_state.abs.sph.cpp ${GEOAC_DIR}geoac.params.sph.cpp ${GEOAC_DIR}geoac.eigenray.sph.mpi.cpp ${GEOAC_DIR}geoac.eqset.sph.strat.cpp ${UTILS_DIR}rk4solver.cpp ${UTILS_DIR}waveforms.cpp ${UTILS_DIR}multitasking.cpp ${GEOAC_DIR}geoac.interface.sph.cpp ${MAINS_DIR}main.infraga.sph.mpi.cpp ${CFLAGS} -o bin/infraga-accel-sph

infraga-accel-sph-rngdep:
	${mpiCC}  ${ATMO_DIR}atmo_io.sph.rngdep.mpi.cpp ${ATMO_DIR}atmo_state.sph.rngdep.cpp ${UTILS_DIR}globe.cpp ${UTILS_DIR}fileIO.cpp ${UTILS_DIR}interpolation.cpp ${ATMO_DIR}atmo_state.abs.sph.cpp ${GEOAC_DIR}geoac.params.sph.cpp ${GEOAC_DIR}geoac.eigenray.sph.mpi.cpp ${GEOAC_DIR}geoac.eqset.sph.rngdep.cpp ${UTILS_DIR}rk4solver.cpp ${GEOAC_DIR}geoac.interface.sph.cpp ${UTILS_DIR}waveforms.cpp ${UTILS_DIR}multitasking.cpp ${MAINS_DIR}main.infraga.sph.rngdep.mpi.cpp ${CFLAGS} -o bin/infraga-accel-sph-rngdep

# Clean the built binaries

clean: 
	rm bin/infraga-2d bin/infraga-3d bin/infraga-3d-rngdep bin/infraga-sph bin/infraga-sph-rngdep

clean_accel: 
	rm bin/infraga-accel-3d bin/infraga-accel-3d-rngdep bin/infraga-accel-sph bin/infraga-accel-sph-rngdep

# Install the executables into the INSTALL_DIR
install: all
	install -m 0755 bin/infraga-2d ${INSTALL_DIR}
	install -m 0755 bin/infraga-3d ${INSTALL_DIR}
	install -m 0755 bin/infraga-3d-rngdep ${INSTALL_DIR}
	install -m 0755 bin/infraga-sph ${INSTALL_DIR}
	install -m 0755 bin/infraga-sph-rngdep ${INSTALL_DIR}	

install_accel: accel
	install -m 0755 bin/infraga-accel-3d ${INSTALL_DIR}
	install -m 0755 bin/infraga-accel-3d-rngdep ${INSTALL_DIR}
	install -m 0755 bin/infraga-accel-sph ${INSTALL_DIR}
	install -m 0755 bin/infraga-accel-sph-rngdep ${INSTALL_DIR}

# Uninstall the methods from INSTALL_DIRS

uninstall:
	rm ${INSTALL_DIR}/infraga-2d 
	rm ${INSTALL_DIR}/infraga-3d 
	rm ${INSTALL_DIR}/infraga-3d-rngdep 
	rm ${INSTALL_DIR}/infraga-sph
	rm ${INSTALL_DIR}/infraga-sph-rngdep 

uninstall_accel:
	rm ${INSTALL_DIR}/infraga.3d.accel
	rm ${INSTALL_DIR}/infraga.3d.rngdep.accel
	rm ${INSTALL_DIR}/infraga.sph.accel
	rm ${INSTALL_DIR}/infraga.sph.rngdep.accel

