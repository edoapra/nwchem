# $Id$

  OBJ = cneo_bas.o 			cneo_constrain.o \
		cneo_dft_cvg.o 		cneo_dft_fockbld.o \
		cneo_dft_info.o 	 \
		cneo_dft.o 			cneo_dftmain0.o \
		cneo_fock_file.o \
		cneo_fock.o 		cneo_get_info.o \
		cneo_input.o \
		cneo_schwarz.o \
		cneo_tddft_energy.o \
		cneo_tddft_init.o   cneo_tdgrad_compute_g.o \
		cneo_tdgrad_compute_w.o 	cneo_tdgrad_init.o \
		cneo_tdgrad_product.o 	    cneo_tdgrad_response.o \
		cneo_tdgrad_solve_init.o    cneo_tdgrad_solve.o \
		cneo_tdgrad_tidy.o          cneo_tdgrad.o \
		get_j_dm.o 			util_newton_solve.o

  OBJ_OPTIMIZE = \
		cneo_diis_driver.o  \
		cneo_grad.o   cneo_dft_scf.o

      LIBRARY =	libcneo.a

    USES_BLAS = get_j_dm.F   cneo_fock_file.F   cneo_tdgrad_compute_w.F   cneo_tdgrad_compute_g.F   cneo_grad.F

LIB_DEFINES = 

LIB_INCLUDES = -I../nwdft/include -I../nwdft/grid -I../ddscf -I../nwdft/lr_tddft_grad

HEADERS = ccneo.fh ccneo_constrain.fh ccneo_dft.fh ccneo_schwarz.fh ccneo_tdgrad_solve.fh

include ../config/makefile.h
ifeq ($(_USE_SCALAPACK),Y)
_USE_ELPA := $(shell ${GA_PATH}/bin/ga-config  --use_elpa| awk ' /1/ {print "Y"}')

ifeq ($(_USE_ELPA),Y)
 LIB_DEFINES += -DUSE_ELPA
endif
#assuming we have GAPZHEEVD
GOTGAPZHEEVD= Y
ifeq ($(GOTGAPZHEEVD),Y)
 LIB_DEFINES += -DUSE_GAPZHEEVD
endif
_PDSYEVR := $(shell [ -e ${NWCHEM_TOP}/src/tools/install/lib/libga.a ] &&  nm ${NWCHEM_TOP}/src/tools/install/lib/libga.a|egrep T |egrep  ga_pdsyevr |awk '/ga_pdsyevr/ {print "Y"}')
ifeq ($(_PDSYEVR),Y) 
 LIB_DEFINES += -DUSE_PDSYEVR
endif
endif
include ../config/makelib.h