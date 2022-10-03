#!/bin/bash

# get the tool which aids making wheel package
# it will automatically generate python wrapper for decorated (doxygen) C code
rm -rf cBinder || exit 1
if [ -d "$DIR" ]; then
    echo "cBinder already exists"
    cd cBinder || exit 1
    git pull || exit 1
else
    git clone https://github.com/Tetrite/cBinder.git || exit 1
fi

# copy source files and headers to local folder, so the docker container can access them
rm -rf libamtrack || exit 1
mkdir libamtrack || exit 1
cp -rf ../../../.git libamtrack/  || exit 1 # needed to evaluate version number from GIT tags
cp -rf ../../../include libamtrack/ || exit 1
cp -rf ../../../src libamtrack/ || exit 1

# get necessary dependency (latest) version of GSL library
if [ -d "gsl-latest" ]; then
    echo "gsl-latest already exists"
else
    rm -f gsl.tar.gz || exit 1
    wget --quiet "http://ftpmirror.gnu.org/gnu/gsl/gsl-latest.tar.gz" --output-document=gsl.tar.gz || exit 1

    # unpack into gsl-x.y directory
    rm -rf gsl-?.* || exit 1
    tar -zxf gsl.tar.gz || exit 1

    # move into gsl-latest directory
    rm -rf gsl-latest || exit 1
    mv gsl-?.* gsl-latest || exit 1
fi


# make directory for GSL installation
if [ -d "gsl" ]; then
    echo "gsl already exists"
else
    mkdir gsl || exit 1
fi

# run the docker process which using manylinux1 machine (SLC6 based) will compile
# the C extensions using old enough libraries
# in such way the will be as portable as possible
rm -rf generated || exit 1
docker run --rm -e PLAT=manylinux1_x86_64 -v `pwd`:/io quay.io/pypa/manylinux1_x86_64 /io/generate_inside_docker.sh || exit 1

# list the generated files
ls -alh generated/dist/
