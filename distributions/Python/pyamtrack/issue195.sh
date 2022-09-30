#!/bin/bash
sudo ./make_wheel_package.sh
python -m pip install --force-reinstall  generated/dist/wheelhouse/pyamtrack-0.13.0-py3-none-manylinux_2_5_x86_64.manylinux1_x86_64.whl
python issue195.py