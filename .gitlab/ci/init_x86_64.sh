#!/bin/bash

source ./build_flavor.sh

# Determine is you have CVMFS installed
if [ ! -d "/cvmfs" ]; then
    echo "No CVMFS detected, please install it."
    exit 1
fi

if [ ! -d "/cvmfs/clicdp.cern.ch" ]; then
    echo "No clicdp CVMFS repository detected, please add it."
    exit 1
fi

if [ ! -d "/cvmfs/sft.cern.ch" ]; then
    echo "No sft CVMFS repository detected, please add it."
    exit 1
fi

# General variables
CLICREPO=/cvmfs/clicdp.cern.ch
SFTREPO=/cvmfs/sft.cern.ch


#--------------------------------------------------------------------------------
#     Compiler
#--------------------------------------------------------------------------------

if [ ${COMPILER_TYPE} == "gcc" ]; then
    source ${CLICREPO}/compilers/gcc/7.1.0/x86_64-${OS}/setup.sh
fi
if [ ${COMPILER_TYPE} == "llvm" ]; then
    source ${CLICREPO}/compilers/llvm/4.0.0/x86_64-${OS}/setup.sh
fi

#--------------------------------------------------------------------------------
#     Python
#--------------------------------------------------------------------------------
export PYTHONDIR=${CLICREPO}/software/Python/2.7.13/${BUILD_FLAVOUR}
export PATH=${PYTHONDIR}/bin:$PATH
export LD_LIBRARY_PATH=${PYTHONDIR}/lib:${LD_LIBRARY_PATH}

#--------------------------------------------------------------------------------
#     CMake
#--------------------------------------------------------------------------------

export CMAKE_HOME=${CLICREPO}/software/CMake/3.8.1/${BUILD_FLAVOUR}
export PATH=${CMAKE_HOME}/bin:$PATH

#--------------------------------------------------------------------------------
#     ROOT
#--------------------------------------------------------------------------------

export ROOTSYS=${CLICREPO}/software/ROOT/6.08.06/${BUILD_FLAVOUR}
export PYTHONPATH="$ROOTSYS/lib:$PYTHONPATH"
export PATH="$ROOTSYS/bin:$PATH"
export LD_LIBRARY_PATH="$ROOTSYS/lib:$LD_LIBRARY_PATH"
export CMAKE_PREFIX_PATH="$ROOTSYS:$CMAKE_PREFIX_PATH"

#--------------------------------------------------------------------------------
#     Geant4
#--------------------------------------------------------------------------------

export G4INSTALL=${CLICREPO}/software/Geant4/10.03.p01/${BUILD_FLAVOUR}
export G4LIB=$G4INSTALL/lib64/Geant4-10.3.1/
export G4ENV_INIT="${G4INSTALL}/bin/geant4.sh"
export G4SYSTEM="Linux-g++"
export CMAKE_PREFIX_PATH="$G4INSTALL:$CMAKE_PREFIX_PATH"

#--------------------------------------------------------------------------------
#     Ninja
#--------------------------------------------------------------------------------

export Ninja_HOME=${CLICREPO}/software/Ninja/1.7.2/${BUILD_FLAVOUR}
export PATH="$Ninja_HOME:$PATH"

#--------------------------------------------------------------------------------
#     Eigen
#--------------------------------------------------------------------------------

export Eigen_HOME=${CLICREPO}/software/Eigen/3.3.3/${BUILD_FLAVOUR}/
export Eigen3_DIR=${Eigen_HOME}/share/eigen3/cmake/
export CMAKE_PREFIX_PATH="$Eigen3_DIR:$CMAKE_PREFIX_PATH"

#--------------------------------------------------------------------------------
#     Doxygen
#--------------------------------------------------------------------------------

export Doxygen_HOME=${SFTREPO}/lcg/releases/doxygen/1.8.11-ae1d3/${BUILD_FLAVOUR}/bin/
export PATH="$Doxygen_HOME:$PATH"

#--------------------------------------------------------------------------------
#     Git
#--------------------------------------------------------------------------------

export Git_HOME=${CLICREPO}/software/git/2.13.2/${BUILD_FLAVOUR}
export PATH=${Git_HOME}/bin:${PATH}

#--------------------------------------------------------------------------------
#     LCIO
#--------------------------------------------------------------------------------

export LCIO=${CLICREPO}/software/LCIO/2.8.0/${BUILD_FLAVOUR}/
export CMAKE_PREFIX_PATH="$LCIO:$CMAKE_PREFIX_PATH"
