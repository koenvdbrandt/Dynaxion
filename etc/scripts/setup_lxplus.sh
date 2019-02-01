if [ ! -d "/cvmfs/sft.cern.ch" ]; then
    echo "CVMFS not available"
    return
fi

if [ "$(uname)" == "Linux" ]; then
    if [ "$( cat /etc/*-release | grep Scientific )" ]; then
        OS=slc6
    elif [ "$( cat /etc/*-release | grep CentOS )" ]; then
        OS=centos7
    else
        echo "Cannot detect OS, falling back to SLC6"
        OS=slc6
    fi
else
    echo "Unknown OS"
    exit 1
fi

# Get our directory and load the CI init
ABSOLUTE_PATH=`dirname $(readlink -f ${BASH_SOURCE[0]})`

# Set the compiler type to LLVM to also load clang-format and clang-tidy
export COMPILER_TYPE="llvm"

# Load default configuration
source $ABSOLUTE_PATH/../../.gitlab/ci/init_x86_64.sh

# Load different Geant for the moment, because CLICdp version does not have QT
# FIXME: This should not be a fixed directory

if [ ${OS} == "slc6"]; then
    export G4INSTALL=/cvmfs/sft.cern.ch/lcg/releases/Geant4/10.04.p02-fd180/x86_64-slc6-gcc7-opt/
    export CMAKE_PREFIX_PATH="/cvmfs/sft.cern.ch/lcg/releases/clhep/2.4.1.0-2c56f/x86_64-slc6-gcc7-opt:$CMAKE_PREFIX_PATH"
elif [ ${OS} == "centos7"]; then
    export G4INSTALL=/cvmfs/sft.cern.ch/lcg/releases/Geant4/10.04.p02-fe3eb/x86_64-centos7-gcc7-opt/
    export CMAKE_PREFIX_PATH="/cvmfs/sft.cern.ch/lcg/releases/clhep/2.4.1.0-2c56f/x86_64-centos7-gcc7-opt:$CMAKE_PREFIX_PATH"
else
    echo "Could not find SFT Geant4 version for ${OS}"
    exit 1
fi

source $G4INSTALL/bin/geant4.sh
export CMAKE_PREFIX_PATH="$G4INSTALL:$CMAKE_PREFIX_PATH"
# Point to corresponding CLHEP
