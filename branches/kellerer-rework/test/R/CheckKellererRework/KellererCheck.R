rm(list=ls())
library(lattice)

df		<-	read.table("KellererCheck.csv", header = T, sep = ";")
df$value	<-	df$lowest.left.limit * df$step^df$bin.no

xyplot(	log10(frequency) ~ bin.no, #log10(value),
		df,
		type 	= 'p',
		groups = n.conv,
		auto.key = list(space = "right"))




