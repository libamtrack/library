N2	<-	5
step	<-	2^(1/N2)
U	<-	log(step)


LH	<-	10
LF	<-	7

FLF	<-	1/U*log(exp(U*LH)-exp(U*LF))

FLF	<-	log(step^LH - step^LF)/log(step)

FLF	<-	LH + log(1.0 - exp(-U*(LH - LF)))/U

FLF	<-	LH + log(1 - step^(LF-LH)) / log(step)