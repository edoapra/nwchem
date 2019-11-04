#!/bin/bash
TARBALL=https://github.com/nwchemgit/nwchem/releases/download/v7.0.0-beta1/nwchem-7.0.0-release.revision-5bcf0416-src.2019-11-01.tar.bz2
if [[ ! -z "$TARBALL" ]]; then
    cd $TRAVIS_BUILD_DIR/..
    pwd
    mv nwchem nwchem.git
    curl -LJ $TARBALL -o nwchem.tar.bz2
    tar xjf nwchem.tar.bz2
    ln -sf nwchem-7.0.0 nwchem
fi
# source env. variables
 source $TRAVIS_BUILD_DIR/travis/nwchem.bashrc
 ls -lrt $TRAVIS_BUILD_DIR|tail -2
 cd $TRAVIS_BUILD_DIR/src
     make nwchem_config
if [[ -z "$TARBALL" ]]; then
 if [[ "$USE_64TO32" == "y" ]]; then
     echo " CONVERSION 64_to_32"
     make 64_to_32 >& 6log &
 fi
fi
 env
