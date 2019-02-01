#!/bin/bash

# Determine which OS you are using
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

# Determine which compiler to use
if [ -z ${COMPILER_TYPE} ]; then
    echo "No compiler type set, falling back to GCC."
    COMPILER_TYPE="gcc"
fi
if [ ${COMPILER_TYPE} == "gcc" ]; then
    echo "Compiler type set to GCC."
    COMPILER_VERSION="gcc7"
fi
if [ ${COMPILER_TYPE} == "llvm" ]; then
    echo "Compiler type set to LLVM."
    COMPILER_VERSION="llvm40"
fi

# Choose build type
if [ -z ${BUILD_TYPE} ]; then
    BUILD_TYPE=opt
fi

export BUILD_FLAVOUR=x86_64-${OS}-${COMPILER_VERSION}-${BUILD_TYPE}
