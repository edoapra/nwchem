#!/usr/bin/env bash

if [[  -z "${NWCHEM_TOP}" ]]; then
    dir3=$(dirname `pwd`)
    dir2=$(dirname "$dir3")
    NWCHEM_TOP=$(dirname "$dir2")
fi

if [[ "$BLAS_SIZE" == 4 ]] && [[ -z "$USE_64TO32"   ]] ; then
    if [[ "$NWCHEM_TARGET" != "LINUX" ]] && [[ "$NWCHEM_TARGET" != "MACX" ]] ; then
    echo USE_64TO32 must be set when BLAS_SIZE=4 on 64-bit architectures
    exit 1
    fi
fi
if [[ ! -z "$BUILD_OPENBLAS"   ]] ; then
    BLASOPT="-L`pwd`/../lib -lnwc_openblas"
fi
plumed_branch=cvhd
github_user=edoapra
git clone --depth=1 -b $plumed_branch https://github.com/$github_user/plumed2 plumed2
cd plumed2
if [[  -z "${FC}" ]]; then
    FC=gfortran
fi
FC_EXTRA=$(${NWCHEM_TOP}/src/config/strip_compiler.sh ${FC})
LDFLAGS_EXTRA="-L/usr/local/opt/fftw/lib "
if [[  "${FC_EXTRA}" == "gfortran" ]]; then
    LDFLAGS_EXTRA+=" -L"`gfortran -print-file-name=libgfortran.a|sed -e s/libgfortran.a//`" -lgfortran"
fi
if [[ "$BLAS_SIZE" == 8 ]];  then
   ILP64="--enable-ilp64 "
else
    ILP64=" "
fi
echo LDFLAGS_EXTRA is "$LDFLAGS_EXTRA"
./configure --enable-modules=+cvhd --disable-mpi \
	    "$ILP64" \
	    LDFLAGS="$BLASOPT $LDFLAGS_EXTRA" \
	    --prefix=$NWCHEM_TOP/src/libext
make  -j4 
if [[ "$?" != "0" ]]; then
    echo " "
    echo "compilation failed"
    echo " "
    exit 1
fi
make install
cd ..
#cp lib/libscalapack.a ../../../lib/libnwc_scalapack.a
