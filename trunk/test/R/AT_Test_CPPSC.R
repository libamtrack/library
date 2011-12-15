################################################################################################
# R test script for CPPSC
################################################################################################
# Copyright 2006, 2010 The libamtrack team
# 
# This file is part of the AmTrack program (libamtrack.sourceforge.net).
#
# AmTrack is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# AmTrack is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# long with AmTrack (file: copying.txt).
# If not, see <http://www.gnu.org/licenses/>
#
################################################################################################

# clear workspace
rm( list = ls() )

# Load local version of library
require(libamtrack, lib.loc = ".")

require(lattice)

####################### CPPSC function test #########################################

# Mixed field of 3 particles:

E.MeV.u                <- c( 50)
particle.no            <- c( 6012)
fluence.cm2.or.dose.Gy <- c( -10 )                                       # in Gy
material.no            <- 1                                                    # Water, Liquid
RDD.model              <- 3                                                    # Geiss
RDD.parameters         <- 5e-8                                                 # a0
ER.model               <- 4                                                    # Geiss
gamma.model            <- 2                                                    # General hit/target
gamma.parameters       <- c(1,10,1,1,0)                                        # Exp.-sat. (one-hit/one-target) with 10 Gy sat.-dose
N2                     <- 10                                                   # 10 bins / factor of 2
source.no              <- 0                                                    # PSTAR

CPPSC.res              <- AT.run.CPPSC.method(        E.MeV.u                     = E.MeV.u,
									particle.no                 = particle.no,
									fluence.cm2.or.dose.Gy      = fluence.cm2.or.dose.Gy,
									material.no                 = material.no,
									stopping.power.source.no    = source.no,
									rdd.model                   = RDD.model,
									rdd.parameters              = RDD.parameters,
									er.model                    = ER.model,
									gamma.model                 = gamma.model,
									gamma.parameters            = gamma.parameters,
									N2                          = N2,
									fluence.factor              = 1.0,
									write.output                = F,
									shrink.tails                = T,
									shrink.tails.under          = 1e-30,
									adjust.N2                   = T,
									lethal.events.mode          = F)

print(CPPSC.res)

