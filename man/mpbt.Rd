\name{mpbt}
\alias{mpbt}
\title{Testing and correcting for publication bias of multivariate meta-analysis}
\description{
Testing and correcting for publication bias of multivariate meta-analysis}
\usage{
mpbt(data, method, type, k)
}
\arguments{
  \item{data}{dataset}
  \item{method}{"nn.cl" indicating the score test for detecting publication bias of MMA}
  \item{type}{either "continuous" or "binary" indicating the type of outcomes}
  \item{k}{integer indicating the number of outcomes}
}
\details{
This function returns the test statistics for testing publication bias of multivariate
meta-analysis using regression method.} 
\section{A score test for detecting publication bias in multivariate meta-analysis}{
Publication bias occurs when the publication of research results depends not only on the quality of the research but also on the direction, magnitude or statistical significance of the results. 
The consequence is that published studies may not represent all valid studies undertaken, 
and this bias may threaten the validity of systematic reviews and meta-analyses - on which comparative effectiveness research and evidence-based medicine increasingly relies. 
Multivariate meta-analysis has recently received increasing attention for its potential ability in reducing bias and improving statistical efficiency by borrowing information across outcomes. 
However, detecting and accounting for publication bias are more challenging in multivariate meta-analysis setting because some studies may be completely unpublished whereas other studies may selectively report part of multiple outcomes. 
Hong et al. (2015) propose a pseudolikelihood-based score test for detecting publication bias in multivariate random-effects meta-analysis. 
This is the first test for detecting publication bias in multivariate meta-analysis setting.
}

\value{
\code{mpbt.TS} returns the test statistic and p value of the score test.
}
\references{Hong, C., Chu, H. and Chen Y. (2015). A score test for detecting publication bias in multivariate random-effects meta-analysis (in preparation).
}
\author{Chuan Hong}

\examples{
data(prostate)
fit.mpbt=mpbt(data=prostate, method = "nn.cl", type = "continuous", k=2)
summary(fit.mpbt)
}
\keyword{Multivariate meta-analysis}
\keyword{Publication bias}
\keyword{Score test}
