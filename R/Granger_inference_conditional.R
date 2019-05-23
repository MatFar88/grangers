#' Inference on conditional Granger-causality
#'
#' \verb{Granger.inference.conditional} provides bootstrap inference for the Granger-causality
#' conditional spectrum of a time series \verb{x} (effect variable) on a time series \verb{z} (conditioning variable)
#' respect to a time series \verb{y} (cause variable). It requires packages \href{https://CRAN.R-project.org/package=vars}{vars}  
#' and \href{https://CRAN.R-project.org/package=tseries}{tseries}.
#' @param x univariate time series.
#' @param y univariate time series (of the same length of \verb{x}).
#' @param z univariate time series (of the same length of \verb{x}).
#' @param ic.chosen estimation method parameter \verb{ic}
#' to be passed to function \link[vars]{VAR} of package \href{https://CRAN.R-project.org/package=vars}{vars}.
#' Defaults to ''SC'' (Schwarz criterion). Alternatives are \verb{c(''AIC'',''HQ'',''SC'',''FPE'')}.
#' @param max.lag maximum number of lags \verb{lag.max} 
#' to be passed to function \code{\link[vars]{VAR}}.
#' Defaults to \verb{min(4, length(x) - 1)}.
#' @param nboots number of bootstrap series to be computed by function \code{\link[tseries]{tsbootstrap}}
#' 	 of package \href{https://CRAN.R-project.org/package=tseries}{tseries}. It defaults to 1000.
#' @param conf prescribed confidence level. It defaults to 0.95.
#' @param bp matrix containing previously simulated bootstrap series, having as rows time points, as columns variables \verb{x} and \verb{y} (in this order). It defaults to NULL.
#' @param p1 parameter \verb{p} to be passed to function \link[vars]{VAR}.
#'	It corresponds to the number of lags of the first VAR model. Defaults to 0.
#' @param p2 parameter \verb{p} to be passed to function \link[vars]{VAR}.
#'	It corresponds to the number of lags of the second VAR model. Defaults to 0.
#' @param ts_boot boolean equal to 1 if the stationary bootstrap of 
#' 	Politis and Romano (1994) is applied, 0 otherwise. It defaults to 1.
#' @description Inference on the conditional Granger-causality spectrum is provided generating
#' 	bootstrap time series by the stationary boostrap of Politis and Romano (1994).
#' 	For computational details we refer to Ding et al. (2006) and Farne' and Montanari (2018).
#' @return \verb{frequency}: frequencies used by Fast Fourier Transform.
#' @return \verb{n}: time series length.
#' @return \verb{nboots}: number of bootstrap series used.
#' @return \verb{confidence_level}: prescribed confidence level.
#' @return \verb{stat_yes}: boolean equal to 0 if no stationary VAR 
#'		is estimated across bootstrap samples, 1 otherwise.
#' @return \verb{non_stationarity_rate_1}: percentage of non-stationary VAR models (at
#' 	least one root larger than one) estimated on bootstrapped \verb{x} and \verb{z}.
#' @return \verb{non_stationarity_rate_2}: percentage of non-stationary VAR models (at
#' 	least one root larger than one) estimated on bootstrapped \verb{x} and \verb{y} and \verb{z}.
#' @return \verb{delay1_mean}: mean number of delays of stationary VAR models estimated on \verb{x} and \verb{z}.
#' @return \verb{delay2_mean}: mean number of delays of stationary VAR models estimated on \verb{x} and \verb{y} and \verb{z}.
#' @return \verb{quantile_conditional_causality_y.to.x.on.z}: computed quantile of the Granger-
#' 	causality conditional spectrum from \verb{y} to \verb{x} on \verb{z}. Differently from function
#' 	\code{\link[grangers]{Granger.inference.unconditional}}, this function provides only the quantile
#' 	of the unidirectional causality from \verb{y} to \verb{x}.
#' @return \verb{freq_y.to.x.on.z}: frequencies at which the Granger-causality conditional spectrum
#'	from \verb{y} to \verb{x} condtional on \verb{z} exceeds the computed threshold.
#' @return \verb{q_max_x.on.z}: computed quantile of the Granger-
#' 	causality conditional spectrum from \verb{y} to \verb{x} on \verb{z} under Bonferroni correction. Differently from function
#' 	\code{\link[grangers]{Granger.inference.unconditional}}, this function provides only the quantile
#' 	of the unidirectional causality from \verb{y} to \verb{x}.
#' @return \verb{freq_max_y.to.x.on.z}: frequencies at which the Granger-causality conditional spectrum
#'	from \verb{y} to \verb{x} conditional on \verb{z} exceeds the computed threshold under Bonferroni correction.
#' @return The result is returned invisibly if plot is TRUE.
#' @author Matteo Farne', Angela Montanari, \email{matteo.farne2@@unibo.it}
#' @seealso \link[vars]{VAR} and \code{\link[tseries]{tsbootstrap}}.
#' @examples
#' Granger.inference.conditional(RealGdp.rate.ts,m3.rate.ts,hicp.rate.ts,plot=T)
#' bp_all<-array(0,dim=c(length(RealGdp.rate.ts),nboots,3))
#' bp_all[,,1]<-tsbootstrap(RealGdp.rate.ts,nb=nboots)
#' bp_all[,,2]<-tsbootstrap(m3.rate.ts,nb=nboots)
#' bp_all[,,3]<-tsbootstrap(hicp.rate.ts,nb=nboots)
#' bp_used<-bp_all[,,c(1,2,3)]
#' inf_cond_m3.to.gdp.by.hicp_0.95<-Granger.inference.conditional(RealGdp.rate.ts,m3.rate.ts,hicp.rate.ts,bp=bp_used)
#' inf_cond_m3.to.gdp.by.hicp_0.90<-
#' Granger.inference.conditional(RealGdp.rate.ts,m3.rate.ts,hicp.rate.ts,conf=0.90,bp=bp_used)
#' @references Politis D. N. and Romano  J. P., (1994). ''The Stationary
#'    Bootstrap''. \emph{Journal of the American Statistical Association}, 89, 1303--1313.
#' @references Ding, M., Chen, Y., Bressler, S.L., 2006. Granger Causality: Basic Theory and
#' 	Application to Neuroscience, Chap.17. \emph{Handbook of Time Series Analysis
#' 	Recent Theoretical Developments and Applications}.
#' @references Farne', M., Montanari, A., 2018. A bootstrap test to detect prominent Granger-causalities across frequencies. 
#'	\emph{Submitted}.
#' @export

Granger.inference.conditional<-function (x, y, z, ic.chosen = "SC", max.lag = min(4, length(x) - 
    1), plot = F, type.chosen = "none", p1=0, p2=0, nboots = 1000, conf = 0.95, 
    bp = NULL,ts_boot=1) 
{

    if (length(x) == 1) {
        return("The length of x is only 1")
    }
    if (length(x) != length(y)) {
        return("x and y do not have the same length")
    }
    if (length(x) != length(z)) {
        return("x and z do not have the same length")
    }
    if (max.lag > length(x) - 1) {
        return("The chosen number of lags is larger than or equal to the time length")
    }

	if (p1==0){
	model1=VAR(cbind(x,z),ic=ic.chosen,lag.max=max.lag,type=type.chosen)
	}

	if (p1>0){
	model1=VAR(cbind(x,z),p=p,type=type.chosen)
	}
	
	if (p2==0){
	model2=VAR(cbind(x,y,z),ic=ic.chosen,lag.max=max.lag,type=type.chosen)
	}

	if (p2>0){
	model2=VAR(cbind(x,y,z),p=p,type=type.chosen)
	}

    if(ts_boot==1){
    freq.good = spec.pgram(y, plot = F)$freq/frequency(y)
    if (is.array(bp) != TRUE) {
        if (!("tseries" %in% installed.packages())) {
            install.packages("tseries")
            library(tseries)
        }
        x_bp <- tsbootstrap(x, nb = nboots)
        y_bp <- tsbootstrap(y, nb = nboots)
        z_bp <- tsbootstrap(z, nb = nboots)
    }
    if (is.array(bp) == TRUE) {
        x_bp <- bp[, , 1]
        y_bp <- bp[, , 2]
        z_bp <- bp[, , 3]

    }
    freq.good = spec.pgram(y, plot = F)$freq/frequency(y)
    }

      delay1_bp <- vector("numeric", nboots)
      delay2_bp <- vector("numeric", nboots)

    test_stationarity_1 <- vector("numeric", nboots)
    test_stationarity_2 <- vector("numeric", nboots)
    top_bp_y.to.x.on.z <- vector("numeric", nboots)
    freq.curr.l<-vector("numeric", nboots)
    stat_rate <- vector("numeric", nboots)
    cause_bp_y.to.x.on.z <- array(0, dim = c(nboots, length(freq.good)))

    for (w in 1:nboots) {
  if(ts_boot==1){
  xz_mat<-as.data.frame(cbind(x_bp[, w], z_bp[, w]));
  xyz_mat<-as.data.frame(cbind(x_bp[, w],y_bp[, w], z_bp[, w]));
  colnames(xz_mat)<-c("x_bp","z_bp")
  colnames(xyz_mat)<-c("x2_bp","y2_bp","z2_bp")
  if(p1>0 && p2>0){
  model1_bp<-VAR(xz_mat, type=type.chosen,p=model1$p)
  model2_bp<-VAR(xyz_mat, type=type.chosen,p=model2$p)
  G.xy <- Granger.conditional(xyz_mat[, 1], xyz_mat[, 2], xyz_mat[, 3], plot=F, type=type.chosen, p1=model1$p,p2=model2$p)
  }
  if(p1==0 && p2==0){
  model1_bp<-VAR(xz_mat, ic=ic.chosen, 
            lag.max=max.lag, type=type.chosen)
  model2_bp<-VAR(xyz_mat, ic=ic.chosen, 
            lag.max=max.lag, type=type.chosen)
  G.xy <- Granger.conditional(xyz_mat[, 1], xyz_mat[, 2], xyz_mat[, 3], ic.chosen, 
            max.lag, F, type.chosen)
  }
  }
    ##
  delay1_bp[w]<-model1_bp$p   
  delay2_bp[w]<-model2_bp$p  
  ##
	  freq.curr.l[w]<-length(G.xy$Conditional_causality_y.to.x.on.z)
        cause_bp_y.to.x.on.z[w,] <- G.xy$Conditional_causality_y.to.x.on.z
        if (length(which(abs(G.xy$roots_1) >= 1)) > 0) {
            test_stationarity_1[w] = 1
        }
        if (length(which(abs(G.xy$roots_2) >= 1)) > 0) {
            test_stationarity_2[w] = 1
        }
    }
    cause_bp_y.to.x.on.z<-cause_bp_y.to.x.on.z[,1:min(freq.curr.l)]
    for (w in 1:(nboots)) {
        	top_bp_y.to.x.on.z[w] <- median(cause_bp_y.to.x.on.z[w, ])
	
    }

 
    stationary <- intersect(which(test_stationarity_1 == 0), 
        which(test_stationarity_2 == 0))
    q_x.on.z <- quantile(top_bp_y.to.x.on.z[stationary], conf)
    alpha_bonf<-(1-conf)/length(freq.good)
    conf_bonf<-1-alpha_bonf
    q_max_x.on.z <- quantile(top_bp_y.to.x.on.z[stationary], conf_bonf)
    non_stationarity_rate_1 <- sum(test_stationarity_1)/nboots
    non_stationarity_rate_2 <- sum(test_stationarity_2)/nboots
    stat_rate<-length(stationary)/nboots


    if (length(stationary)>=nboots/nboots){
    stat_yes=1;
    
    n <- G.xy$n

    GG_x.on.z<-Granger.conditional(x,y,z,ic.chosen,max.lag,F)$Conditional_causality_y.to.x.on.z
    signif_x.on.z<-which(GG_x.on.z>q_x.on.z)
    signif_max_x.on.z<-which(GG_x.on.z>q_max_x.on.z)

    GG <- list(freq.good, n, nboots, conf, stat_yes,  non_stationarity_rate_1, 
        non_stationarity_rate_2, mean(delay1_bp[test_stationarity_1==0]),mean(delay2_bp[test_stationarity_2==0]), q_x.on.z, freq.good[signif_x.on.z],q_max_x.on.z,freq.good[signif_max_x.on.z])
    names(GG) <- c("frequency", "n", "nboots","confidence_level","stat_yes", "non_stationarity_rate_1", "non_stationarity_rate_2", "delay1_mean","delay2_mean","quantile_conditional_causality_y.to.x.on.z", "freq_y.to.x.on.z","q_max_x.on.z","freq_max_y.to.x.on.z")
    }

    if (length(stationary)<nboots/nboots){
    stat_yes=0;
    GG<-list(stat_yes)
    names(GG)<-c("stat_yes")
    }
	
    if (plot == F) {
        return(GG)
    }
    if (plot == T) {
        par(mfrow = c(1, 1))
        plot(freq.good, GG_x.on.z, type = "l", main = "Conditional Granger-causality y to x on z")
        abline(h = q_x.on.z)
	  abline(h = q_max_x.on.z)
    }
}