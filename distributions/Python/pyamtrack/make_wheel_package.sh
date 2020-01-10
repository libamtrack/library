#!/bin/bash

mkdir workspace
cp run_me_inside_docker.sh workspace
docker run --rm -e PLAT=manylinux1_x86_64 -v workspace:/io quay.io/pypa/manylinux1_x86_64 run_me_inside_docker.sh

