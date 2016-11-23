
filename="full_callset.tcga.22.fscl.gpu"

df <- read.table(filename, header=FALSE)
y <- sort(df$V7)
n <- length(y)
x <- sort(-log10((1:n)/n))
m <- ceiling(max(x))
plot(x,y,pch=19,xlim=c(0,m),ylim=c(0,m))
lines(c(0,m),c(0,m),col="red")

error.quantile <- 0.95
qm <- 10^m
alpha <- qm:1
beta <- qm - alpha + 1
xtop <- qbeta(error.quantile, alpha, beta)
xbot <- qbeta(1.0 - error.quantile, alpha, beta)
lines(-log10(alpha/(qm+1)),-log10(xbot), lty=2)
lines(-log10(alpha/(qm+1)),-log10(xtop), lty=2)
