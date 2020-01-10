#!/bin/bash

#mkdir workspace
#cp run_me_inside_docker.sh workspace
git clone https://github.com/Tetrite/cBinder.git

mkdir libamtrack
cp -r ../../../src libamtrack/
cp -r ../../../include libamtrack/

wget http://ftp.task.gda.pl/pub/gnu/gsl/gsl-latest.tar.gz
tar -zxvf gsl-latest.tar.gz

docker run --rm -e PLAT=manylinux1_x86_64 -v `pwd`:/io quay.io/pypa/manylinux1_x86_64 /io/run_me_inside_docker.sh

ls -alh generated/dist/