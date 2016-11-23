

filename="full_callset.tcga.22.fscl"
output.filename=""

args <- commandArgs(trailingOnly=TRUE)
i <- 1
while(i <= length(args)) {
    if (args[i] == '-f') {
        filename = args[i+1];
        i <- i+1
    } else if (args[i] == '-o') {
        output.filename = args[i+1];
        i <- i+1
    } else {
        cat(sprintf("Command line option %s not recognized\n",args[i]))
    }
    i <- i+1
}

out.pdf <- 0
if (output.filename != "") {
    out.pdf <- 1
    cairo_pdf(file=output.filename, width=8.5, height=11, onefile=TRUE, pointsize=12, antialias="gray")
}


df <- read.table(filename, header=FALSE)
y <- sort(df$V7)
n <- length(y)
x <- sort(-log10((1:n)/n))
m <- ceiling(max(x))
plot(x,y,pch=19,xlim=c(0,m),ylim=c(0,m), xlab="-log10(p) expected quantiles", ylab="-log10(p) observed quantiles")
lines(c(0,m),c(0,m),col="red")

error.quantile <- 0.95
qm <- 10^m
alpha <- qm:1
beta <- qm - alpha + 1
xtop <- qbeta(error.quantile, alpha, beta)
xbot <- qbeta(1.0 - error.quantile, alpha, beta)
lines(-log10(alpha/(qm+1)),-log10(xbot), lty=2)
lines(-log10(alpha/(qm+1)),-log10(xtop), lty=2)

if (out.pdf) {
    dev.off()
}
