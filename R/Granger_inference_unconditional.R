#' Inference on unconditional Granger-causality
#'
#' \verb{Granger.inference.unconditional} provides bootstrap inference for the Granger-causality
#' 	unconditional spectrum of a time series \verb{x} (effect variable) respect to a time series
#' 	\verb{y} (cause variable). It requires packages \href{https://CRAN.R-project.org/package=vars}{vars} and \href{https://CRAN.R-project.org/package=tseries}{tseries}.
#' @param x univariate time series.
#' @param y  univariate time series (of the same length of \verb{x}).
#' @param ic.chosen estimation method parameter \verb{ic} to be passed to function \link[vars]{VAR} of
#' 	package \href{https://CRAN.R-project.org/package=vars}{vars}. Defaults to ''SC'' (Schwarz criterion). Alternatives are \verb{c(''AIC'',''HQ'',''SC'',''FPE'')}.
#' @param max.lag maximum number of lags \verb{lag.max} to be passed to function \code{\link[vars]{VAR}}.
#' 	Defaults to \verb{min(4, length(x) - 1)}.
#' @param plot logical; if TRUE, it returns the plot of unconditional Granger-causality
#' 	spectra on both directions with computed thresholds. Defaults to FALSE.
#' @param type.chosen parameter \verb{type} to be passed to function \code{\link[vars]{VAR}}.
#'	Defaults to \verb{''none''}. Alternatives are \verb{c(''none'',''const'',''trend'')}.
#' @param p parameter \verb{p} to be passed to function \link[vars]{VAR}.
#'	Defaults to 0.
#' @param nboots number of bootstrap series to be computed by function \code{\link[tseries]{tsbootstrap}}
#' 	 of package \href{https://CRAN.R-project.org/package=tseries}{tseries}. It defaults to 1000.
#' @param conf prescribed confidence level. It defaults to 0.95.
#' @param bp matrix containing previously simulated bootstrap series, having as rows
#' 	time points, as columns variables \verb{x} and \verb{y} (in this order). It defaults to NULL.
#' @param ts_boot boolean equal to 1 if the stationary bootstrap of 
#' 	Politis and Romano (1994) is applied, 0 otherwise. It defaults to 1.
#' @description Inference on the unconditional Granger-causality spectrum is provided generating
#' 	bootstrap time series by the stationary boostrap of Politis and Romano (1994).
#' 	For computational details we refer to Ding et al. (2006) and Farne' and Montanari (2018).
#' @return \verb{frequency}: frequencies used by Fast Fourier Transform.
#' @return \verb{n}: time series length.
#' @return \verb{nboots}: number of bootstrap series used.
#' @return \verb{confidence_level}: prescribed confidence level.
#' @return \verb{stat_yes}: boolean equal to 0 if no stationary VAR 
#'		is estimated across bootstrap samples, 1 otherwise.
#' @return \verb{non_stationarity_rate}: percentage of non-stationary VAR models (at
#' 	least one root larger than one) estimated on bootstrapped \verb{x} and \verb{y}.
#' @return \verb{delay_mean}: mean number of delays of stationary VAR models estimated on \verb{x} and \verb{y}.
#' @return \verb{quantile_unconditional_causality_y.to.x}: computed quantile of the Granger-causality
#'    	unconditional spectrum from \verb{y} to \verb{x}.
#' @return \verb{quantile_unconditional_causality_x.to.y}: computed quantile of the Granger-causality
#'    	unconditional spectrum from \verb{x} to \verb{y}.
#' @return \verb{freq_y.to.x}: frequencies at which the Granger-causality unconditional spectrum
#'		from \verb{y} to \verb{x} exceeds the computed threshold.
#' @return \verb{freq_x.to.y}: frequencies at which the Granger-causality unconditional spectrum
#'		from \verb{x} to \verb{y} exceeds the computed threshold.
#' @return \verb{q_max_x}: computed quantile of the Granger-causality
#'    	unconditional spectrum from \verb{y} to \verb{x} under Bonferroni correction.
#' @return \verb{q_max_y}: computed quantile of the Granger-causality
#'    	unconditional spectrum from \verb{x} to \verb{y} under Bonferroni correction.
#' @return \verb{freq_max_y.to.x}: frequencies at which the Granger-causality unconditional spectrum
#'		from \verb{y} to \verb{x} exceeds the computed threshold under Bonferroni correction.
#' @return \verb{freq_max_x.to.y}: frequencies at which the Granger-causality unconditional spectrum
#'		from \verb{x} to \verb{y} exceeds the computed threshold under Bonferroni correction.
#' @return The result is returned invisibly if plot is TRUE.
#' @author Matteo Farne', Angela Montanari, \email{matteo.farne2@@unibo.it}
#' @seealso \link[vars]{VAR} and \code{\link[tseries]{tsbootstrap}}.
#' @examples
#' 	RealGdp.rate.ts<-euro_area_indicators[,1]
#'	m3.rate.ts<-euro_area_indicators[,2]
#' 	inf_uncond_m3_0.95<-Granger.inference.unconditional(RealGdp.rate.ts,m3.rate.ts,nboots=10)
#' @references Politis D. N. and Romano  J. P., (1994). ''The Stationary
#'    Bootstrap''. \emph{Journal of the American Statistical Association}, 89, 1303--1313.
#' @references Ding, M., Chen, Y., Bressler, S.L., 2006. Granger Causality: Basic Theory and
#' 	Application to Neuroscience, Chap.17. \emph{Handbook of Time Series Analysis
#' 	Recent Theoretical Developments and Applications}.
#' @references Farne', M., Montanari, A., 2018. A bootstrap test to detect prominent Granger-causalities across frequencies. 
#'	<arXiv:1803.00374>, \emph{Submitted}.
#' @export
#' @import vars tseries
#' @importFrom graphics abline par
#' @importFrom stats coef frequency median pf qf quantile residuals spec.pgram
#' @importFrom utils install.packages installed.packages

Granger.inference.unconditional<-function (x, y, ic.chosen = "SC", max.lag = min(4, length(x) - 
    1), plot = F, type.chosen = "none",p=0,nboots = 1000, conf = 0.95, 
    bp = NULL,ts_boot=1) 
{

    if (length(x) == 1) {
        return("The length of x is only 1")
    }
    if (length(x) != length(y)) {
        return("x and y do not have the same length")
    }
    if (max.lag > length(x) - 1) {
        return("The chosen number of lags is larger than or equal to the time length")
    }

if(!require("vars")){
return("The packages 'vars' could not be found. Please install it to 
proceed.")
}

if(!require("tseries")){
return("The packages 'tseries' could not be found. Please install it to 
proceed.")
}

require(vars)
require(tseries)

    if(ts_boot==1){
    if (is.array(bp) != TRUE) {

        x_bp <- tsbootstrap(x, nb = nboots)
        y_bp <- tsbootstrap(y, nb = nboots)
    }
    if (is.array(bp) == TRUE) {
        x_bp <- bp[, , 1]
        y_bp <- bp[, , 2]
    }
    freq.good = spec.pgram(y, plot = F)$freq/frequency(x)
    }

    if (p==0){
       mod=VAR(cbind(x,y),ic=ic.chosen,lag.max=max.lag,type.chosen)
    }

    if (p>0){
       mod=VAR(cbind(x,y),ic=ic.chosen,lag.max=max.lag,type.chosen,p=p)
    }

    delay_bp <- vector("numeric", nboots)
    test_stationarity <- vector("numeric", nboots)
    top_bp_y.to.x <- vector("numeric", nboots)
    top_bp_x.to.y <- vector("numeric", nboots)
    cause_y.to.x_bp <- array(0, dim = c(nboots, length(freq.good)))
    cause_x.to.y_bp <- array(0, dim = c(nboots, length(freq.good)))
    for (w in 1:nboots) {
  xy_mat<-as.data.frame(cbind(x_bp[, w], y_bp[, w]));
  colnames(xy_mat)<-c("x_bp","y_bp")
  if(p>0){
  mod_bp<-VAR(xy_mat, ic=ic.chosen, lag.max=max.lag, type.chosen,p=mod$p)
  G.xy <- Granger.unconditional(xy_mat[, 1], xy_mat[, 2], plot=F, type.chosen, p=mod$p)
  }
  if(p==0){
  mod_bp<-VAR(xy_mat, ic=ic.chosen, 
            lag.max=max.lag, type.chosen)
  G.xy <- Granger.unconditional(xy_mat[, 1], xy_mat[, 2], ic.chosen, 
            max.lag, F, type.chosen)
  }
  delay_bp[w]<-mod_bp$p        
  cause_y.to.x_bp[w, ] <- G.xy$Unconditional_causality_y.to.x
  cause_x.to.y_bp[w, ] <- G.xy$Unconditional_causality_x.to.y
        if (length(which(abs(G.xy$roots) >= 1)) > 0) {
            test_stationarity[w] = 1
        }
    }

    for (w in 1:(nboots)) {
        top_bp_y.to.x[w] <- median(cause_y.to.x_bp[w, ])
        top_bp_x.to.y[w] <- median(cause_x.to.y_bp[w, ])
    }

    q_x <- quantile(top_bp_y.to.x[test_stationarity == 0], conf)
    q_y <- quantile(top_bp_x.to.y[test_stationarity == 0], conf)

    alpha_bonf<-(1-conf)/length(freq.good)
    conf_bonf<-1-alpha_bonf

    q_max_x <- quantile(top_bp_y.to.x[test_stationarity == 0], conf_bonf)
    q_max_y <- quantile(top_bp_x.to.y[test_stationarity == 0], conf_bonf)

    non_stationarity_rate <- sum(test_stationarity)/nboots

    stat_rate<-1-non_stationarity_rate

    if (stat_rate>=1/nboots){
    stat_yes=1;
    n <- G.xy$n    
    GG_x_orig<-Granger.unconditional(x,y,ic.chosen,max.lag,F)$Unconditional_causality_y.to.x
    GG_y_orig<-Granger.unconditional(x,y,ic.chosen,max.lag,F)$Unconditional_causality_x.to.y
    signif_x<-which(GG_x_orig>q_x)
    signif_y<-which(GG_y_orig>q_y)
    signif_max_x<-which(GG_x_orig>q_max_x)
    signif_max_y<-which(GG_y_orig>q_max_y)
    GG <- list(freq.good, n, nboots, conf, stat_yes, non_stationarity_rate, mean(delay_bp[test_stationarity == 0]), q_x, q_y, 
        freq.good[signif_x],freq.good[signif_y],q_max_x,q_max_y, freq.good[signif_max_x], freq.good[signif_max_y])
    names(GG) <- c("frequency", "n", "nboots", "confidence_level", "stat_yes",  "non_stationarity_rate", "delay_mean", "quantile_unconditional_causality_y.to.x", 
        "quantile_unconditional_causality_x.to.y",
	  "freq_y.to.x","freq_x.to.y","q_max_x","q_max_y","freq_max_y.to.x","freq_max_x.to.y")
    
    }

    if (stat_rate==0){
    stat_yes=0;
    GG<-list(stat_yes)
    names(GG)<-c("stat_yes")
    }
    
    if (plot == F) {
        return(GG)
    }
    if (plot == T) {
        par(mfrow = c(1, 2))
        plot(freq.good, GG_x_orig, type = "l", main = "Unconditional Granger-causality y to x")
        abline(h = q_x)
        plot(freq.good, GG_y_orig, type = "l", main = "Unconditional Granger-causality x to y")
        abline(h = q_y)
    }
}