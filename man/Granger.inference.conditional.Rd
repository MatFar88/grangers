% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Granger_inference_conditional.R
\name{Granger.inference.conditional}
\alias{Granger.inference.conditional}
\title{Inference on conditional Granger-causality}
\usage{
Granger.inference.conditional(x, y, z, ic.chosen = "SC",
  max.lag = min(4, length(x) - 1), plot = F, type.chosen = "none",
  p1 = 0, p2 = 0, nboots = 1000, conf = 0.95, bp = NULL,
  ts_boot = 1)
}
\arguments{
\item{x}{univariate time series.}

\item{y}{univariate time series (of the same length of \verb{x}).}

\item{z}{univariate time series (of the same length of \verb{x}).}

\item{ic.chosen}{estimation method parameter \verb{ic}
to be passed to function \link[vars]{VAR} of package \href{https://CRAN.R-project.org/package=vars}{vars}.
Defaults to ''SC'' (Schwarz criterion). Alternatives are \verb{c(''AIC'',''HQ'',''SC'',''FPE'')}.}

\item{max.lag}{maximum number of lags \verb{lag.max} 
to be passed to function \code{\link[vars]{VAR}}.
Defaults to \verb{min(4, length(x) - 1)}.}

\item{plot}{logical; if TRUE, it returns the plot of unconditional Granger-causality
spectra on both directions with computed thresholds. Defaults to FALSE.}

\item{type.chosen}{parameter \verb{type} to be passed to function \code{\link[vars]{VAR}}.
Defaults to \verb{''none''}. Alternatives are \verb{c(''none'',''const'',''trend'')}.}

\item{p1}{parameter \verb{p} to be passed to function \link[vars]{VAR}.
It corresponds to the number of lags of the first VAR model. Defaults to 0.}

\item{p2}{parameter \verb{p} to be passed to function \link[vars]{VAR}.
It corresponds to the number of lags of the second VAR model. Defaults to 0.}

\item{nboots}{number of bootstrap series to be computed by function \code{\link[tseries]{tsbootstrap}}
of package \href{https://CRAN.R-project.org/package=tseries}{tseries}. It defaults to 1000.}

\item{conf}{prescribed confidence level. It defaults to 0.95.}

\item{bp}{matrix containing previously simulated bootstrap series, having as rows time points, as columns variables \verb{x} and \verb{y} (in this order). It defaults to NULL.}

\item{ts_boot}{boolean equal to 1 if the stationary bootstrap of 
Politis and Romano (1994) is applied, 0 otherwise. It defaults to 1.}
}
\value{
\verb{frequency}: frequencies used by Fast Fourier Transform.

\verb{n}: time series length.

\verb{nboots}: number of bootstrap series used.

\verb{confidence_level}: prescribed confidence level.

\verb{stat_yes}: boolean equal to 0 if no stationary VAR 
	is estimated across bootstrap samples, 1 otherwise.

\verb{non_stationarity_rate_1}: percentage of non-stationary VAR models (at
	least one root larger than one) estimated on bootstrapped \verb{x} and \verb{z}.

\verb{non_stationarity_rate_2}: percentage of non-stationary VAR models (at
	least one root larger than one) estimated on bootstrapped \verb{x} and \verb{y} and \verb{z}.

\verb{delay1_mean}: mean number of delays of stationary VAR models estimated on \verb{x} and \verb{z}.

\verb{delay2_mean}: mean number of delays of stationary VAR models estimated on \verb{x} and \verb{y} and \verb{z}.

\verb{quantile_conditional_causality_y.to.x.on.z}: computed quantile of the Granger-
	causality conditional spectrum from \verb{y} to \verb{x} on \verb{z}. Differently from function
	\code{\link[grangers]{Granger.inference.unconditional}}, this function provides only the quantile
	of the unidirectional causality from \verb{y} to \verb{x}.

\verb{freq_y.to.x.on.z}: frequencies at which the Granger-causality conditional spectrum
from \verb{y} to \verb{x} condtional on \verb{z} exceeds the computed threshold.

\verb{q_max_x.on.z}: computed quantile of the Granger-
	causality conditional spectrum from \verb{y} to \verb{x} on \verb{z} under Bonferroni correction. Differently from function
	\code{\link[grangers]{Granger.inference.unconditional}}, this function provides only the quantile
	of the unidirectional causality from \verb{y} to \verb{x}.

\verb{freq_max_y.to.x.on.z}: frequencies at which the Granger-causality conditional spectrum
from \verb{y} to \verb{x} conditional on \verb{z} exceeds the computed threshold under Bonferroni correction.

The result is returned invisibly if plot is TRUE.
}
\description{
Inference on the conditional Granger-causality spectrum is provided generating
	bootstrap time series by the stationary boostrap of Politis and Romano (1994).
	For computational details we refer to Ding et al. (2006) and Farne' and Montanari (2018).
}
\details{
\verb{Granger.inference.conditional} provides bootstrap inference for the Granger-causality
conditional spectrum of a time series \verb{x} (effect variable) on a time series \verb{z} (conditioning variable)
respect to a time series \verb{y} (cause variable). It requires packages \href{https://CRAN.R-project.org/package=vars}{vars}  
and \href{https://CRAN.R-project.org/package=tseries}{tseries}.
}
\examples{
	RealGdp.rate.ts<-euro_area_indicators[,1]
m3.rate.ts<-euro_area_indicators[,2]
hicp.rate.ts<-euro_area_indicators[,4]
	inf_cond_m3.to.gdp.by.hicp_0.95<-
Granger.inference.conditional(RealGdp.rate.ts,m3.rate.ts,hicp.rate.ts,nboots=10)
}
\references{
Politis D. N. and Romano  J. P., (1994). ''The Stationary
   Bootstrap''. \emph{Journal of the American Statistical Association}, 89, 1303--1313.

Ding, M., Chen, Y., Bressler, S.L., 2006. Granger Causality: Basic Theory and
	Application to Neuroscience, Chap.17. \emph{Handbook of Time Series Analysis
	Recent Theoretical Developments and Applications}.

Farne', M., Montanari, A., 2018. A bootstrap test to detect prominent Granger-causalities across frequencies. 
<arXiv:1803.00374>, \emph{Submitted}.
}
\seealso{
\link[vars]{VAR} and \code{\link[tseries]{tsbootstrap}}.
}
\author{
Matteo Farne', Angela Montanari, \email{matteo.farne2@unibo.it}
}
