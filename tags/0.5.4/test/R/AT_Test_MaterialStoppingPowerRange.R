################################################################################################
# R test script for material properties, stopping power, and range functions
################################################################################################
# Copyright 2006, 2011 The libamtrack team
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

# Materials

material.nos       <- c(1,2,3,4,5,6,7,11,12)
particle.nos       <- AT.particle.no.from.particle.name("1H")

E.MeV.u            <- 10^seq(from = 0, to = 3, length.out = 40)

df                 <- expand.grid( E.MeV.u                        = E.MeV.u,
                                   particle.no                    = particle.nos,
                                   material.no                    = material.nos,
                                   Stopping.Power.PSTAR.MeV.g.cm2 = 0,
                                   Stopping.Power.Bethe.MeV.g.cm2 = 0,
                                   Stopping.Power.Ratio.Bethe     = 0)

#df$material.name   <- AT.material.name.from.material.no(material.no = df$material.no)


# Get Bethe Stopping Power
for (i in 1:nrow(df)){
    # DEBUG i <- 2
    df$Stopping.Power.Bethe.MeV.g.cm2[i]   <- AT.Stopping.Power.Mass.Bethe.MeV.cm2.g(         E.MeV.u          = df$E.MeV.u[i],
                                                                                              particle.no      = df$particle.no[i],
                                                                                              material.no      = df$material.no[i],
                                                                                              E.restricted.keV = 0)[1]
}

for( cur.mat.no in unique(df$material.no)){
    ii                                 <-  df$material.no == cur.mat.no
    jj                                 <-  df$material.no == 1
    df$Stopping.Power.Ratio.Bethe[ii]  <-  df$Stopping.Power.Bethe.MeV.g.cm2[jj] / df$Stopping.Power.Bethe.MeV.g.cm2[ii]
}


xyplot(   log10(Stopping.Power.Bethe.MeV.g.cm2) ~ log10(E.MeV.u),
          df,
          type         = 'l',
          groups       = material.no,
          auto.key     = list(space = 'right'))


xyplot(   df$Stopping.Power.Ratio.Bethe ~ log10(E.MeV.u),
          df,
          type         = 'l',
          groups       = material.no,
          auto.key     = list(space = 'right'))

