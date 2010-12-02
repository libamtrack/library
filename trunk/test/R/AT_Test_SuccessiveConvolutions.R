# Testing script for AT_SuccessiveConvolutions routine
# Created: 2010-10-05
# Creator: greilich

rm(list = ls())

# load libAmTrack library
try(dyn.load("../../wrapper/R/R_direct_access/libamtrack.dll"))
try(dyn.load("../../wrapper/R/R_direct_access/libamtrack.so"))
try(dyn.load("../../wrapper/R/R_direct_access/libamtrack.dylib"))

# load wrapping scripts
source("../../wrapper/R/R_direct_access/libamtrack.R")

library(lattice)

u					<- 10
n.bins.f			<- 100
N2					<- 10
n.bins.f.used		<- 10

df					<- data.frame(f.d.Gy = 2^seq(1, 10, length.out = 100))
df$f.dd.Gy			<- c(2.0 - 1.888749, diff(df$f.d.Gy))
df$f				<- c(rep(0.05,10), rep(0, 90))

f0					<- 0.5

write.output		<- F
shrink.tails		<- T
shrink.tails.under 	<- 1e-30
adjust.N2			<- T
		
res	<- AT.SC.SuccessiveConvolutions(	u		= u,
								n.bins.f		= n.bins.f,
								N2				= N2,
								n.bins.f.used	= n.bins.f.used,
								f.d.Gy			= df$f.d.Gy,
								f.dd.Gy			= df$f.dd.Gy,
								f				= df$f,
								write.output	= write.output,
								shrink.tails	= shrink.tails,
								shrink.tails.under = shrink.tails.under,
								adjust.N2		= adjust.N2)
								

df$which <- "start"
res$f$which <- "final"

df.plot <- data.frame(	f.d.Gy	= c(df$f.d.Gy, res$f$f.d.Gy),
						f		= c(df$f, res$f$f),
						which	= c(df$which, res$f$which))

p1 <- xyplot(f ~ log10(f.d.Gy),
		df.plot,
		type = 's',
		groups = which,
		auto.key = T)

						
pdf("AT_test_successive_convolutions.pdf")						

p1

dev.off()