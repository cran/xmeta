\name{summary.galaxy}
\alias{summary.galaxy}
\title{Summarize the objects \code{galaxy}}
\description{
    Summary a model of class \code{galaxy} fitted by \code{galaxy}.
}
\usage{
      \method{summary}{galaxy}(object,...) 
}

\arguments{  
    \item{object}{an object inheriting from class \code{galaxy}.}
	 \item{...}{ additional arguments; currently none is used.}
}


\value{
  A list with the following components: adjusted parameter estimates, variances of estimated parameters
}

\references{
Chen, Y., Hong, C., Chu, H., (2015). Galaxy plot and a multivariate method for correcting publication bias in multivariate meta-analysis (in preparation).

}



\seealso{\code{\link{galaxy}}}
\examples{
data(prostate)
fit.galaxy=galaxy(data=prostate, type = "continuous",
 method="galaxy.cl", k=2, L=1, estimator="R0", maxiter=150)
summary(fit.galaxy)
}


\keyword{summary}
