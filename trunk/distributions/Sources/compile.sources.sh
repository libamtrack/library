#!/bin/sh

autoreconf --force --install
./configure
make all
make dist

