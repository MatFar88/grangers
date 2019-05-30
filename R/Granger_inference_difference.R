#' Inference on the difference between unconditional and conditional Granger-causality
#'
#' \verb{Granger.inference.difference} provides bootstrap inference for the difference between
#' 	the Granger-causality unconditional spectrum of a time series \verb{x} (effect variable)
#' 	respect to a time series \verb{y} (cause variable) and the Granger-causality conditional
#' 	spectrum of a time series \verb{x} (effect variable) on a time series \verb{z} (conditioning variable)
#' 	respect to a time series \verb{y} (cause variable). It requires packages \href{https://CRAN.R-project.org/package=vars}{vars} and \href{https://CRAN.R-project.org/package=tseries}{tseries}.
#' @param x univariate time series.
#' @param y  univariate time series (of the same length of \verb{x}).
#' @param z  univariate time series (of the same length of \verb{x}).
#' @param ic.chosen estimation method parameter \verb{ic} to be passed to function \link[vars]{VAR} of
#' 	package ''vars''. Defaults to ''SC'' (Schwarz criterion). Alternatives are \verb{c(''AIC'',''HQ'',''SC'',''FPE'')}.
#' @param max.lag maximum number of lags \verb{lag.max} to be passed to function \code{\link[vars]{VAR}}.
#' 	Defaults to \verb{min(4, length(x) - 1)}.
#' @param plot logical; if TRUE, it returns the plot of the difference between the unconditional
#' 	Granger-causality spectrum from \verb{y} to \verb{x} and the conditional Granger-causality
#' 	spectrum from \verb{y} to \verb{x} on \verb{z} with upper and lower computed thresholds.
#'	Defaults to FALSE.
#' @param type.chosen parameter \verb{type} to be passed to function \code{\link[vars]{VAR}}.
#'	Defaults to \verb{''none''}. Alternatives are \verb{c(''none'',''const'',''trend'')}.
#' @param p parameter \verb{p} to be passed to function \link[vars]{VAR}.
#'	It corresponds to the number of delays for unconditional GC. Defaults to 0.
#' @param p1 parameter \verb{p} to be passed to function \link[vars]{VAR}.
#'	It corresponds to the number of lags of the first VAR model. Defaults to 0.
#' @param p2 parameter \verb{p} to be passed to function \link[vars]{VAR}.
#' @param nboots number of bootstrap series to be computed by function \code{\link[tseries]{tsbootstrap}}
#' 	 of package \href{https://CRAN.R-project.org/package=tseries}{tseries}. It defaults to 1000.
#' @param conf prescribed confidence level. It defaults to 0.95.
#' @param bp_orig matrix containing previously simulated bootstrap series, having as rows
#' 	time points, as columns variables \verb{x} and \verb{y} (in this order). It defaults to NULL.
#' @param ts_boot boolean equal to 1 if the stationary bootstrap of 
#' 	Politis and Romano (1994) is applied, 0 otherwise. It defaults to 1.
#' @description Inference on the difference between unconditional and conditional Granger-causality
#' 	spectrum is provided generating bootstrap time series by the stationary boostrap of
#' 	Politis and Romano (1994).
#' 	For computational details we refer to Ding et al. (2006) and Farne' and Montanari (2018).
#' @return \verb{frequency}: frequencies used by Fast Fourier Transform.
#' @return \verb{n}: time series length.
#' @return \verb{nboots}: number of bootstrap series used.
#' @return \verb{confidence_level}: prescribed confidence level.
#' @return \verb{stat_yes}: boolean equal to 0 if no stationary VAR 
#'		is estimated across bootstrap samples, 1 otherwise.
#' @return \verb{non_stationarity_rate}: percentage of estimated non-stationary VAR models (at
#' 	least one root larger than one) on bootstrapped \verb{x} and {y}.
#' @return \verb{non_stationarity_rate_1}: percentage of estimated non-stationary VAR models (at
#' 	least one root larger than one) on bootstrapped \verb{x} and {z}.
#' @return \verb{non_stationarity_rate_2}: percentage of estimated non-stationary VAR models (at
#' 	least one root larger than one) on bootstrapped \verb{x}, \verb{y} and {z}.
#' @return \verb{quantile_difference_inf}: lower computed quantile of the difference between the
#' 	Granger-causality unconditional spectrum from \verb{y} to \verb{x} and the Granger-causality
#' 	conditional spectrum from \verb{y} to \verb{x} on \verb{z}.
#' @return \verb{quantile_difference_sup}: upper computed quantile of the difference between the
#' 	Granger-causality unconditional spectrum from \verb{y} to \verb{x} and the Granger-causality
#' 	conditional spectrum from \verb{y} to \verb{x} on \verb{z}.
#' @return \verb{freq_inf}: frequencies at which the difference between the Granger-causality unconditional spectrum 
#'	from \verb{y} to \verb{x} and the Granger-causality conditional spectrum
#' 	from \verb{y} to \verb{x} on \verb{z} exceeds the lower computed threshold. 
#' @return \verb{freq_sup}: frequencies at which the difference between the Granger-causality unconditional spectrum 
#'	from \verb{y} to \verb{x} and the Granger-causality conditional spectrum
#' 	from \verb{y} to \verb{x} on \verb{z} exceeds the upper computed threshold. 
#' @return \verb{quantile_difference_max_inf}: lower computed quantile of the difference between the
#' 	Granger-causality unconditional spectrum from \verb{y} to \verb{x} and the Granger-causality
#' 	conditional spectrum from \verb{y} to \verb{x} on \verb{z} under Bonferroni correction.
#' @return \verb{quantile_difference_max_sup}: upper computed quantile of the difference between the
#' 	Granger-causality unconditional spectrum from \verb{y} to \verb{x} and the Granger-causality
#' 	conditional spectrum from \verb{y} to \verb{x} on \verb{z} under Bonferroni correction.
#' @return \verb{freq_max_inf}: frequencies at which the difference between the Granger-causality unconditional 
#'	spectrum from \verb{y} to \verb{x} and the Granger-causality conditional spectrum
#' 	from \verb{y} to \verb{x} on \verb{z} exceeds the lower computed threshold under Bonferroni correction.
#' @return \verb{freq_max_sup}: frequencies at which the difference between the Granger-causality unconditional 
#'	spectrum from \verb{y} to \verb{x} and the Granger-causality conditional spectrum
#' 	from \verb{y} to \verb{x} on \verb{z} exceeds the upper computed threshold under Bonferroni correction.
#' @return The result is returned invisibly if plot is TRUE.
#' @author Matteo Farne', Angela Montanari, \email{matteo.farne2@@unibo.it}
#' @seealso \link[vars]{VAR} and \code{\link[tseries]{tsbootstrap}}.
#' @examples
#' 	RealGdp.rate.ts<-euro_area_indicators[,1]
#'	m3.rate.ts<-euro_area_indicators[,2]
#'	hicp.rate.ts<-euro_area_indicators[,4]
#' 	inf_diff_pre_hicp.to.gdp_0.95<-
#' 	Granger.inference.difference(RealGdp.rate.ts,m3.rate.ts,hicp.rate.ts,nboots=10)
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

Granger.inference.difference<-function (x, y, z, ic.chosen = "SC", max.lag = min(4, length(x) - 
    1), plot = F, type.chosen = "none", p=0,p1=0,p2=0,nboots = 1000, conf=0.95, bp_orig = NULL,ts_boot=1)
{


    if (p==0){
       mod=VAR(cbind(x,y),ic=ic.chosen,lag.max=max.lag,type.chosen)
    }

    if (p>0){
       mod=VAR(cbind(x,y),ic=ic.chosen,lag.max=max.lag,type.chosen,p=p)
    }


    if (length(x) == 1) {
        return("The length of x is only 1")
    }
    if (length(x) != length(y)) {
        return("x and y do not have the same length")
    }
    if (max.lag > length(x) - 1) {
        return("The chosen number of lags is larger than or equal to the time length")
    }

    ##


if(!requireNamespace("vars")){
return("The packages 'vars' could not be found. Please install it to 
proceed.")
}

if(!requireNamespace("tseries")){
return("The packages 'tseries' could not be found. Please install it to 
proceed.")
}

requireNamespace("vars")
requireNamespace("tseries")


	if (p1==0){
	model1=VAR(cbind(x,z),ic=ic.chosen,lag.max=max.lag,type.chosen)
	}

	if (p1>0){
	model1=VAR(cbind(x,z),p=p1,type.chosen)
	}
	
	if (p2==0){
	model2=VAR(cbind(x,y,z),ic=ic.chosen,lag.max=max.lag,type.chosen)
	}

	if (p2>0){
	model2=VAR(cbind(x,y,z),p=p2,type.chosen)
	}

if(ts_boot==1){
    if (is.array(bp_orig) != TRUE) {

        x_bp <- tsbootstrap(x, nb = nboots)
        y_bp <- tsbootstrap(y, nb = nboots)
        z_bp <- tsbootstrap(z, nb = nboots)
    }
    if (is.array(bp_orig) == TRUE) {
        x_bp <- bp_orig[, , 1]
        y_bp <- bp_orig[, , 2]
        z_bp <- bp_orig[, , 3]
    }
    freq.good = spec.pgram(y, plot = F)$freq/frequency(y)

    no_freq=0;
}

    test_stationarity <- vector("numeric", nboots)
    cause_y.to.x.on.z_bp <- array(0, dim = c(nboots, length(freq.good)))
    cause_y.to.x_bp <- array(0, dim = c(nboots, length(freq.good)))
    test_stationarity_1 <- vector("numeric", nboots)
    test_stationarity_2 <- vector("numeric", nboots)
    cause_diff_y.to.x.on.z_bp <- array(0, dim = c(nboots, length(freq.good)))
    cause_diff_y.to.x.on.z_bp_signed <- array(0, dim = c(nboots, length(freq.good)))    
    top_diff_y.to.x.on.z_bp_signed <- vector("numeric", nboots)    

    for (w in 1:nboots) {
  xy_mat<-as.data.frame(cbind(x_bp[, w], y_bp[, w]));
  colnames(xy_mat)<-c("x_bp","y_bp")
  if(p>0){
  mod_bp<-VAR(xy_mat, type.chosen,p=mod$p)
  G.xy <- Granger.unconditional(xy_mat[, 1], xy_mat[, 2], plot=F, type.chosen, p=mod$p)
  }
  if(p==0){
  mod_bp<-VAR(xy_mat, ic=ic.chosen, 
            lag.max=max.lag, type.chosen)
  G.xy <- Granger.unconditional(xy_mat[, 1], xy_mat[, 2], ic.chosen, 
            max.lag, F, type.chosen)

  }

  cause_y.to.x_bp[w, ] <- G.xy$Unconditional_causality_y.to.x

        if (length(which(abs(G.xy$roots) >= 1)) > 0) {
            test_stationarity[w] = 1}

  if(ts_boot==1){
  xz_mat<-as.data.frame(cbind(x_bp[, w], z_bp[, w]));
  xyz_mat<-as.data.frame(cbind(x_bp[, w],y_bp[, w], z_bp[, w]));
  colnames(xz_mat)<-c("x_bp","z_bp")
  colnames(xyz_mat)<-c("x_bp","y_bp","z_bp")
  if(p1>0 && p2>0){
  model1_bp<-VAR(xz_mat, type.chosen,p=model1$p)
  model2_bp<-VAR(xyz_mat, type.chosen,p=model2$p)
  GG.xy <- Granger.conditional(xyz_mat[, 1], xyz_mat[, 2], xyz_mat[, 3], plot=F, type.chosen, p1=model1$p,p2=model2$p)
  }
  if(p1==0 && p2==0){
  model1_bp<-VAR(xz_mat, ic=ic.chosen, 
            lag.max=max.lag, type.chosen)
  model2_bp<-VAR(xyz_mat, ic=ic.chosen, 
            lag.max=max.lag, type.chosen)
  GG.xy <- Granger.conditional(xyz_mat[, 1], xyz_mat[, 2], xyz_mat[, 3], ic.chosen, 
            max.lag, F, type.chosen)
  }
  }

	  #if(is.na(GG.xy)[1]==F){
        cause_y.to.x.on.z_bp[w, ] <- GG.xy$Conditional_causality_y.to.x.on.z
        #}

	  if (length(which(abs(GG.xy$roots_1) >= 1)) > 0) {
            test_stationarity_1[w] = 1
        }
        if (length(which(abs(GG.xy$roots_2) >= 1)) > 0) {
            test_stationarity_2[w] = 1
	  }

	  cause_diff_y.to.x.on.z_bp_signed[w,]<-(cause_y.to.x_bp[w, ] - cause_y.to.x.on.z_bp[w,])
	  top_diff_y.to.x.on.z_bp_signed[w] <- median((cause_diff_y.to.x.on.z_bp_signed[w, ]))

}


    stationary <- intersect(which(test_stationarity == 0), intersect(which(test_stationarity_1 == 
        0), which(test_stationarity_2 == 0)))
    non_stationarity_rate <- sum(test_stationarity)/nboots
    non_stationarity_rate_1 <- sum(test_stationarity_1)/nboots
    non_stationarity_rate_2 <- sum(test_stationarity_2)/nboots

    stat_rate<-length(stationary)/nboots

    if (length(stationary)>=nboots/nboots){
    stat_yes=1;

    n <- G.xy$n


    GG_x_orig <- Granger.unconditional(x, y, ic.chosen, max.lag, 
        F)$Unconditional_causality_y.to.x
    GG_x.on.z <- Granger.conditional(x, y, z, ic.chosen, max.lag, 
        F)$Conditional_causality_y.to.x.on.z


    q_diff_x_inf <- quantile(top_diff_y.to.x.on.z_bp_signed[stationary],(1-conf)/2)
    q_diff_x_sup <- quantile(top_diff_y.to.x.on.z_bp_signed[stationary],1-(1-conf)/2) 
    diff_GG <- (GG_x_orig - GG_x.on.z)
    signif_diff_x_sup <- which(diff_GG > q_diff_x_sup)
    signif_diff_x_inf <- which(diff_GG < q_diff_x_inf)

    conf_bonf<-(1-conf)/length(freq.good)

    q_diff_max_inf <- quantile(top_diff_y.to.x.on.z_bp_signed[stationary],(1-conf_bonf)/2)
    q_diff_max_sup <- quantile(top_diff_y.to.x.on.z_bp_signed[stationary],1-(1-conf_bonf)/2) 
    diff_GG <- (GG_x_orig - GG_x.on.z)
    signif_diff_max_sup <- which(diff_GG > q_diff_max_sup)
    signif_diff_max_inf <- which(diff_GG < q_diff_max_inf)


    GG <- list(freq.good, n, nboots, conf, stat_yes, non_stationarity_rate, non_stationarity_rate_1, non_stationarity_rate_2, 
	q_diff_x_sup, q_diff_x_inf, freq.good[signif_diff_x_sup], freq.good[signif_diff_x_inf],q_diff_max_sup, q_diff_max_inf, freq.good[signif_diff_max_sup],freq.good[signif_diff_max_inf])
    names(GG) <- c("frequency", "n", "nboots", "confidence_level", "stat_yes", "non_stationarity_rate", "non_stationarity_rate_1", "non_stationarity_rate_2", 
        "quantile_difference_sup", "quantile_difference_inf", "freq_sup","freq_inf","quantile_difference_max_sup", "quantile_difference_max_inf", "freq_max_sup","freq_max_inf")

    }

    if (length(stationary)<nboots/nboots){
    stat_yes=0;
    stat_rate=0;
    no_freq=0;
    GG<-list(stat_yes,stat_rate,no_freq)
    names(GG)<-c("stat_yes","stat_rate","no_freq")
    }

    if (plot == F) {
        return(GG)
    }
    if (plot == T) {
        par(mfrow = c(1, 1))
        plot(freq.good, diff_GG, type = "l", main = "Difference Unconditional/Conditional")
        abline(h = q_diff_x_sup)
    }
}