#!/bin/bash

# to install newest python on gitpod type:
# pyenv install 3.10.2 && pyenv local 3.10.2


python3 -m venv venv
source venv/bin/activate

python -m pip install -r generator/requirements.txt

python generator/scripts/header_parser.py --namespace ./NAMESPACE --out-dir ./out/ libamtrack "../../include/*.h"
python generator/scripts/prepare.py --desc-path ./DESCRIPTION --version 0.12.0 --out-dir ./out ../../include ../../src
R --vanilla -e 'setwd("out"); roxygen2::roxygenise(load_code = "source", roclets = "namespace")'
R --vanilla -e 'setwd("out"); pkgbuild::build(binary = TRUE, needs_compilation = TRUE)'

deactivate

# to test run
# sudo R --vanilla -e 'install.packages("libamtrack_0.12.0_R_x86_64-pc-linux-gnu.tar.gz")'
# R --vanilla -e 'library(libamtrack); AT.beta.from.E.single(60)'
