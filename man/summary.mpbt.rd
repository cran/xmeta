\name{summary.mpbt}
\alias{summary.mpbt}
\title{Summarize the objects \code{mpbt}}
\description{
    Summary a model of class \code{mpbt} fitted by \code{mpbt}.
}
\usage{
      \method{summary}{mpbt}(object,...) 
}

\arguments{  
    \item{object}{an object inheriting from class \code{mpbt}.}
	 \item{...}{ additional arguments; currently none is used.}
}


\value{
  A list with the following components: test statistics (mpbt) and p-value.
}

\references{
Hong, C., Chu, H. and Chen Y. (2015). A score test for detecting publication bias in multivariate random-effects meta-analysis (in preparation).

}



\seealso{\code{\link{mpbt}}}
\examples{
data(prostate)
fit.mpbt=mpbt(data=prostate, method = "nn.cl", type = "continuous", k=2)
summary(fit.mpbt)
}


\keyword{summary}
