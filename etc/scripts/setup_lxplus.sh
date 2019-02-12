if [ ! -d "/cvmfs/sft.cern.ch" ]; then
    echo "CVMFS not available"
    return
fi

# Get our directory and load the CI init
ABSOLUTE_PATH=`dirname $(readlink -f ${BASH_SOURCE[0]})`

# Load default configuration
source $ABSOLUTE_PATH/../../.gitlab/ci/build_flavor.sh

# Load different Geant for the moment, because CLICdp version does not have QT
SFTREPO=/cvmfs/sft.cern.ch/lcg/views/LCG_95
source ${SFTREPO}/${BUILD_FLAVOUR}/setup.sh
