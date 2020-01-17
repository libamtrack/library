#!/bin/bash

# get the tool which aids making wheel package
# it will automatically generate python wrapper for decorated (doxygen) C code
rm -rf cBinder
git clone -b fix/81-libamtrack-python3 https://github.com/Tetrite/cBinder.git

# copy source files and headers to local folder, so the docker container can access them
rm -rf libamtrack
mkdir libamtrack
cp -rf ../../../src libamtrack/
cp -rf ../../../include libamtrack/

# get necessary dependency (latest) version of GSL library
wget -q "http://ftpmirror.gnu.org/gnu/gsl/gsl-latest.tar.gz" -O gsl.tar.gz
tar -zxvf gsl.tar.gz
# make directory for GSL installation
rm -rf gsl
mkdir gsl

# run the docker process which using manylinux1 machine (SLC6 based) will compile
# the C extensions using old enough libraries
# in such way the will be as portable as possible
rm -rf generated
docker run --rm -e PLAT=manylinux1_x86_64 -v `pwd`:/io quay.io/pypa/manylinux1_x86_64 /io/generate_inside_docker.sh

# list the generated files
ls -alh generated/dist/