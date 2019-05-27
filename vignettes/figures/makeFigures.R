probs=seq(0,1,0.01)

##truncated normals
require(tmvtnorm)
require(Hmisc)
y0 <- rtmvnorm(1000, mean=1, lower=0)
y1 <- 2*y0 - 0.13*y0^2
y2 <- 0.7*y1 + 2

Ecdf(y1, lwd=3, col="blue")

##truncated normal plot
pdf("cdf-truncated-normal.pdf")
Ecdf(y0, add=F, col="red", lwd=3, main="CDF of Truncated Normal(1,1)",
     ylab="", xlab="x")
segments(x0=0, y0=0.5, x1=quantile(y0,0.5), y1=0.5, lwd=3, lty=2) #horizonal median
segments(x0=quantile(y0,0.5), y0=0.5, x1=quantile(y0,0.5), y1=0, lwd=3, lty=2) #vertical median
segments(x0=0, y0=0.95, x1=quantile(y0,0.95), y1=0.95, lwd=3, lty=2) #horizontal 95 percentile
segments(x0=quantile(y0,0.95), y0=0.95, x1=quantile(y0,0.95), y1=0, lwd=3, lty=2) #vertical 95 percentile
dev.off()
##

##truncated normal quantiles plot
pdf("quantiles-truncated-normal.pdf")
plot(probs, quantile(y0, probs), type="l", col="red", lwd="3",
     main="Quantiles of Truncated Normal(1,1)",
     xlab="tau", ylab="quantiles")
segments(x0=0, y0=quantile(y0,0.5), x1=0.5, y1=quantile(y0,0.5),
         lwd=3, lty=2) #horizontal median
segments(x0=0.5, y0=quantile(y0,0.5), x1=0.5, y1=0,
         lwd=3, lty=2) #vertical median
segments(x0=0, y0=quantile(y0,0.95), x1=0.95, y1=quantile(y0,0.95),
         lwd=3, lty=2) #horizontal 95 percentile
segments(x0=0.95, y0=quantile(y0,0.95), x1=0.95, y1=0,
         lwd=3, lty=2) #vertical 95 percentile
dev.off()
##

##truncated normals qte plot
pdf("cdf-two-distributions.pdf")
Ecdf(y1, add=F, col="blue", lwd=3, main="CDF of Two Distributions",
     ylab="", xlab="x")
Ecdf(y0, add=T, col="red", lwd=3)
segments(x0=quantile(y0,0.5), y0=0.5, x1=quantile(y1,0.5), y1=0.5,
         lwd=3, lty=2) #median QTE
segments(x0=quantile(y0,0.95), y0=0.95, x1=quantile(y1,0.95), y1=0.95,
         lwd=3, lty=2) #95th percentile QTE
dev.off()
##

##
pdf("qte-example.pdf")
plot(probs, quantile(y1, probs)-quantile(y0,probs), type="l", col="purple", lwd="3",
     main="Quantile Treatment Effect",
     xlab="tau", ylab="QTE")
dev.off()
##

## QDiD example
png("qdid-example.png")
yt0 <- rtmvnorm(1000, mean=1, sigma=2)
ecdf.yt0 <- ecdf(yt0)
Ecdf(yt0, xlim=c(-3,6), lwd=3, main="QDiD Counterfactual Distribution",
     ylab="", xlab="")
segments(x0=quantile(yt0,0.5), y0=0.5,
         x1=(quantile(yt0,0.5) + (quantile(y1,0.5) - quantile(y0, 0.5))), y1=0.5,
         lwd=3, lty=2) #median QTE
segments(x0=quantile(yt0,0.95), y0=0.95,
         x1=(quantile(yt0,0.95) + (quantile(y1,0.95) - quantile(y0, 0.95))), y1=0.95,
         lwd=3, lty=2) #median QTE
cf <- yt0 + quantile(y1, ecdf.yt0(yt0)) - quantile(y0, ecdf.yt0(yt0))
Ecdf(cf, add=T, lwd=3, lty=2)
dev.off()
##

## CiC example
##will need to generate two figures for this one
yt0 <- rtmvnorm(1000, mean=1, sigma=2)
ecdf.y0 <- ecdf(y0)
png("cic-example-change.png")
Ecdf(y1, add=F, col="blue", lwd=3, main="CDF of Two Distributions",
     ylab="", xlab="x")
Ecdf(y0, add=T, col="red", lwd=3)
segments(x0=quantile(yt0,0.5), y0=ecdf.y0(quantile(yt0,0.5)),
         x1=quantile(y1,ecdf.y0(quantile(yt0,0.5))),
         y1=ecdf.y0(quantile(yt0,0.5)),
         lwd=3, lty=2) #median QTE
segments(x0=quantile(yt0,0.95), y0=ecdf.y0(quantile(yt0,0.95)),
         x1=quantile(y1,ecdf.y0(quantile(yt0,0.95))),
         y1=ecdf.y0(quantile(yt0,0.95)),
         lwd=3, lty=2) #95th percentile
dev.off()
##second figure
png("cic-example.png")
Ecdf(yt0, xlim=c(-3,6), lwd=3, main="CiC Counterfactual Distribution",
     ylab="", xlab="")
segments(x0=quantile(yt0,0.5), y0=0.5,
         x1=quantile(y1,ecdf.y0(quantile(yt0,0.5))), y1=0.5,
         lwd=3, lty=2) #median QTE
segments(x0=quantile(yt0,0.95), y0=0.95,
         x1=quantile(y1,ecdf.y0(quantile(yt0,0.95))), y1=0.95,
         lwd=3, lty=2) #median QTE
cf <- yt0 + quantile(y1, ecdf.y0(yt0)) - quantile(y0, ecdf.y0(yt0))
Ecdf(cf, add=T, lwd=3, lty=2)
dev.off()
##

##copula picture
png("copula-example.png")
Ecdf(y1, add=F, col="blue", lwd=3, main="Copula of Two Marginal Distributions",
     ylab="", xlab="x")
Ecdf(y0, add=T, col="red", lwd=3)
arrows(x0=quantile(y0,0.7), y0=0.7, x1=quantile(y1,0.7), y1=0.7,
         lwd=3, lty=1, length=0.1) #median QTE
##arrows(x0=quantile(y0,0.76), y0=0.76, x1=quantile(y1,0.98), y1=0.98,
##         lwd=3, lty=1, length=0.1) #
arrows(x0=quantile(y0,0.72), y0=0.72, x1=quantile(y1,0.8), y1=0.8,
         lwd=3, lty=1, length=0.1) #
arrows(x0=quantile(y0,0.74), y0=0.74, x1=quantile(y1,0.65), y1=0.65,
         lwd=3, lty=1, length=0.1) #
##arrows(x0=quantile(y0,0.64), y0=0.64, x1=quantile(y1,0.04), y1=0.04,
##         lwd=3, lty=1, length=0.1) #
arrows(x0=quantile(y0,0.04), y0=0.04, x1=quantile(y1,0.9), y1=0.9,
         lwd=3, lty=1, length=0.1) #
dev.off()
##

##2nd copula picture
png("copula-example-2.png")
Ecdf(y1, add=F, col="blue", lwd=3, main="Copula Stability Assumption",
     ylab="", xlab="x")
Ecdf(y0, add=T, col="red", lwd=3)
Ecdf(y2, add=T, col="green", lwd=3)
arrows(x0=quantile(y1,0.7), y0=0.7, x1=quantile(y2,0.7), y1=0.7,
         lwd=3, lty=1, length=0.1) #median QTE
arrows(x0=quantile(y1,0.72), y0=0.72, x1=quantile(y2,0.8), y1=0.8,
         lwd=3, lty=1, length=0.1) #
arrows(x0=quantile(y1,0.74), y0=0.74, x1=quantile(y2,0.65), y1=0.65,
         lwd=3, lty=1, length=0.1) #
arrows(x0=quantile(y1,0.04), y0=0.04, x1=quantile(y2,0.9), y1=0.9,
         lwd=3, lty=1, length=0.1) #
dev.off()
##

##uniforms
plot(probs, 3*qunif(probs), type="l", lwd=3)
lines(probs, qunif(probs), type="l", lwd=3)

##distribution of uniform
##not right
##plot(probs, 3*punif(probs), type="l", lwd=3)
##lines(probs, punif(probs), type="l", lwd=3)
