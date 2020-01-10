#!/bin/bash

#mkdir workspace
#cp run_me_inside_docker.sh workspace
git clone https://github.com/Tetrite/cBinder.git

docker run --rm -e PLAT=manylinux1_x86_64 -v `pwd`:/io quay.io/pypa/manylinux1_x86_64 /io/run_me_inside_docker.sh

ls -alh workspace
