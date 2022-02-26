#!/bin/bash

python3 -m venv venv
source venv/bin/activate

python -m pip install -r generator/requirements.txt

python prepare.py --desc-path ./DESCRIPTION --version 0.12.0 --out-dir ./out ../../include ../../src

deactivate