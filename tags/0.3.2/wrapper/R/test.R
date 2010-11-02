dyn.load("example.so")
source("example.R")

AT_PrintName()
AT_GetNumber()

n <- 3
E.MeV.u <- c(10.0, 20.0, 100.0)
ranges <- numeric(n)
AT_max_electron_ranges_m( n, E.MeV.u, 1, 3, ranges)

cat( E.MeV.u , "\n")

cat( ranges  , "\n")