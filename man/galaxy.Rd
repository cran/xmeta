\name{galaxy}
\alias{galaxy}
\title{Correcting for publication bias of multivariate meta-analysis}
\description{
Correcting and correcting for publication bias of multivariate meta-analysis}
\usage{
galaxy(data, rhow, type, method, k, L, estimator, side, maxiter)
}
\arguments{
  \item{data}{dataset}
  \item{rhow}{within-study correlation}
  \item{method}{"galaxy.cl" indicating the Galaxy method for detecting and correcting for publication bias of MMA}
  \item{type}{either "continuous" or "binary" indicating the type of outcomes}
  \item{k}{integer indicating the number of outcomes}
  \item{L}{the coefficient for loss function}
  \item{side}{either "left" or "right", indicating on which side of the funnel plot the missing studies should be imputed}
  \item{maxiter}{maximum number of iterations}
  \item{estimator}{either "R0", "L0" or "Q0"}
}
\details{
This function returns the pooled effect size after adjusting for publication bias of multivariate
meta-analysis using the Galaxy method.} 
\value{
\code{res} returns the pooled effect size, covariance matrix, the estimated number of missing and the side of imputation. 
}
\references{Chen, Y., Hong, C., Chu, H., (2015). Galaxy plot and a multivariate method for correcting publication bias in multivariate meta-analysis (in preparation).
}
\author{Yong Chen}


\examples{
library(mvmeta)
y = matrix(0,8,2)
S = inputcov(matrix(runif(8*2,0.1,2),8,2,dimnames=list(NULL,c("V1","V2"))),cor=c(0.5))
Psi = inputcov(1:2,cor=0.3)
mydat=mvmetaSim(y,S,Psi)
y1=mydat[,1]
y2=mydat[,2]
s1=sqrt(S[,1])
s2=sqrt(S[,3])

data=data.frame(y1, s1, y2, s2)
galaxy(data, type = "continuous", method="galaxy.cl", k=2, L=1, estimator="R0", maxiter=150)
}
\keyword{Multivariate meta-analysis}
\keyword{Publication bias}
\keyword{Galaxy method}
