if [ ! -d "/cvmfs/sft.cern.ch" ]; then
    echo "CVMFS not available"
    return
fi

# Get our directory and load the CI init
ABSOLUTE_PATH=`dirname $(readlink -f ${BASH_SOURCE[0]})`

# Load default configuration
source $ABSOLUTE_PATH/../../.gitlab/ci/init_x86_64.sh

# Load different Geant for the moment, because CLICdp version does not have QT
SFTREPO=/cvmfs/sft.cern.ch/lcg/releases

# FIXME: This should not be a fixed directory
if [ ${BUILD_FLAVOUR} == "x86_64-slc6-gcc7-opt" ]; then
    export G4INSTALL=$SFTREPO/Geant4/10.04.p02-0be30/x86_64-slc6-gcc7-opt/
elif [ ${BUILD_FLAVOUR} == "x86_64-centos7-gcc7-opt" ]; then
    export G4INSTALL=$SFTREPO/Geant4/10.04.p02-09680/x86_64-centos7-gcc7-opt/
fi

source $G4INSTALL/bin/geant4.sh
export CMAKE_PREFIX_PATH="$G4INSTALL:$CMAKE_PREFIX_PATH"

# Point to corresponding CLHEP
export CMAKE_PREFIX_PATH="$SFTREPO/clhep/2.4.1.0-2c56f/$BUILD_FLAVOUR:$CMAKE_PREFIX_PATH"
