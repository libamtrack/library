rm(list = ls())
n	<-	10
N2	<-	5
i	<-	n
j	<-	1:n

step	<-	2^(1/N2)

E0	<-	1
F	<-	rep(1.0, n)
D.i	<-	log(1 - step^(j-i))/log(step)
E	<-	E0*step^(n-1)
t	<-	E0*step^(j-1)
E.t	<-	E0*step^(n+D.i-1)
dE	<-	E0*(step^j - step^(j-1))
mid.E	<-	E0*(step^(j+0.5-1))


df	<-	data.frame(n = n, j = j, D.i = D.i, E = E, t = t, E.t = E.t, E.check = t + E.t, mid.E = mid.E, dE = dE, F = F)
df$FDE	<-	df$F * df$dE
df

for (i in 1:n){
	# i <- 3
	sum	<-	0
	for (j in 1:i){
		# j <- 1
		E		<-	E0*step^(i-1)
		t		<-	E0*step^(j-1)
		k		<-	i + D.i[i-j+1]
		E.t		<-	E0*(step^((j-1) - step^(i-1)))
		E.t		<-	E0*step^(i + k-1)
		E.check 	<-	E.t + t
	}
}