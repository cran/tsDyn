% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TVAR.gen.R
\name{TVAR.sim}
\alias{TVAR.sim}
\alias{TVAR.boot}
\title{Simulation of a multivariate Threshold Autoregressive model (TVAR)}
\usage{
TVAR.sim(B, Thresh, nthresh = 1, n = 200, lag = 1, include = c("const",
  "trend", "none", "both"), thDelay = 1, mTh = 1, starting = NULL,
  innov = rmnorm(n, varcov = varcov), varcov = diag(1, nrow(B)),
  show.parMat = FALSE, ...)

TVAR.boot(TVARobject, innov, seed, boot.scheme = c("resample", "wild1",
  "wild2", "check"), ...)
}
\arguments{
\item{B}{Matrix of coefficients to simulate}

\item{Thresh}{The threshold value(s). Vector of length nthresh}

\item{nthresh}{number of threshold (see details)}

\item{n}{Number of observations to create when type="simul"}

\item{lag}{Number of lags to include in each regime}

\item{include}{Type of deterministic regressors to include. NOT WORKING
PROPERLY CURRENTLY if not const}

\item{thDelay}{'time delay' for the threshold variable (as multiple of
embedding time delay d) PLEASE NOTE that the notation is currently different
to univariate models in tsDyn. The left side variable is taken at time t, and
not t+1 as in univariate cases.}

\item{mTh}{combination of variables with same lag order for the transition
variable. Either a single value (indicating which variable to take) or a
combination}

\item{starting}{Starting values (matrix of dimension lag x k). If not given,
set to zero.}

\item{innov}{Innovations used for simulation. Should be matrix of dim n x k.
By default multivariate normal. For the bootstrap case \code{TVAR.boot}, residuals are
resampled if argument is missing.}

\item{varcov}{Variance-covariance matrix for the innovations. By default
identity matrix.}

\item{show.parMat}{Logical. Should the parameter matrix be shown? Useful to
understand how to give right input}

\item{\dots}{Further arguments passed to the underlying (un-exported)
\code{TVAR.gen} function}

\item{TVARobject}{Object of class \code{TVAR} generated by function
\code{\link{TVAR}}}

\item{seed}{Optional. Seed for the random resampling function.}

\item{boot.scheme}{The bootstrap scheme.}
}
\value{
A matrix with the simulated/bootstraped series.
}
\description{
Simulate a multivariate Threshold VAR (TVAR)
}
\details{
This function offers the possibility to generate series following a TVAR.

By giving a matrix of coefficients, on can only simulate a VAR (nthresh=0) or
TVAR (nthresh=1 or 2). One can have a specification with constant (default),
trend, both or none (see arg include). Order in parameters is include/lags
(VECM) and include/lags/include/lags for TVECM, hence, a matrix for a TVECM
with 3 regimes, a const and a 2 lags would have 2 lines and 2*(1+4) columns.
The innovations can be given by the user (a matrix of dim nxk, here n does
not include the starting values!), by default it uses a multivariate normal
distribution, with covariance matrix specified by varcov. The starting values
(of dim lags x k) can be given. The user should take care for their choice,
since it is not sure that the simulated values will cross the threshold even
once.

The matrix \sQuote{B} has to be in the form: constant, trend, lags, then
repeated if many regimes. In case of uncertainty, using
\code{who.parMat=TRUE} will print the matrix as interpreted by the function,
helping the user to feed the right input.

For the bootstrap, the function resamples data from a given TVAR model generated by
\code{\link{TVAR}}, returning the resampled data.
A residual recursive bootstrap is used, where one uses either a simple
resampling, or the Wild bootstrap, either with a normal distribution (wild1) or
inverting the sign randomly (wild2)
}
\examples{


## TVAR.sim: Simulation of a TVAR with 1 threshold
B <- rbind(c(0.11928245, 1.00880447, -0.009974585, -0.089316, 0.95425564, 0.02592617),
        c(0.25283578, 0.09182279,  0.914763741, -0.0530613, 0.02248586, 0.94309347))
colnames(B) <- paste(rep(c("Const", "Lag_1_var1", "Lag_1_var2"), 2), c("Low", "High"), sep="_")
sim <- TVAR.sim(B=B,nthresh=1,n=500, mTh=1, Thresh=5, starting=matrix(c(5.2, 5.5), nrow=1))

#estimate the new serie
TVAR(sim, lag=1, dummyToBothRegimes=TRUE)

## TVAR.boot: Bootstrap a TVAR with two threshold (three regimes)
data(zeroyld)
serie <- zeroyld
mod <- TVAR(data=serie,lag=1, nthresh=1)
TVAR.boot(mod)
}
\seealso{
\code{\link{TVAR}} to estimate the TVAR.  Similar \code{\link{TVECM.sim}} and \code{\link{TVECM.boot}} for \code{\link{TVECM}},
\code{\link{VAR.sim}} and \code{\link{VAR.boot}} for VAR models estimated with \code{\link{lineVar}}.
}
\author{
Matthieu Stigler
}
\keyword{bootstrap}
\keyword{ts}
