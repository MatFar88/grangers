#' Conditional Granger-causality test of Breitung and Candelon (2006)
#'
#' \verb{bc_test_cond} calculates the test of Breitung and Candelon (2006) on the conditional Granger-causality 
#' 	of a time series \verb{x} (effect variable) on a time series \verb{z} (conditioning variable) respect to a time series \verb{y} (cause variable). 
#'	It requires package \href{https://CRAN.R-project.org/package=vars}{vars}.
#' @param x univariate time series.
#' @param y  univariate time series (of the same length of \verb{x}).
#' @param z  univariate time series (of the same length of \verb{x}).
#' @param ic.chosen estimation method parameter \verb{ic} to be passed to function \link[vars]{VAR} of
#' 	package \href{https://CRAN.R-project.org/package=vars}{vars}. Defaults to ''SC'' (Schwarz criterion). Alternatives are \verb{c(''AIC'',''HQ'',''SC'',''FPE'')}.
#' @param max.lag maximum number of lags \verb{lag.max} to be passed to function \code{\link[vars]{VAR}}.
#' 	Defaults to \verb{min(4, length(x) - 1)}.
#' @param plot logical; if TRUE, it returns the plot of conditional Granger-causality
#' 	spectrum. Defaults to FALSE.
#' @param type.chosen parameter \verb{type} to be passed to function \code{\link[vars]{VAR}}.
#' @param p parameter \verb{p} to be passed to function \link[vars]{VAR}.
#'	It corresponds to the number of lags of the second VAR model. Defaults to 0.
#' @param conf prescribed confidence level. It defaults to 0.95.
#' @description Inference on the conditional Granger-causality spectrum is provided by
#' 	the parametric test of Breitung and Candelon (2006).
#' @return \verb{frequency}: frequencies used by Fast Fourier Transform.
#' @return \verb{n}: time series length.
#' @return \verb{confidence_level}: prescribed confidence level.
#' @return \verb{significant_frequencies}: frequencies at which the test is significant..
#' @return \verb{F-test}: computed F-test at each frequency.
#' @return \verb{F-threshold}: F-threshold at each frequency under prescribed confidence level.
#' @return \verb{roots}: roots of the estimated VAR model.
#' @return \verb{delays}: delays of the estimated VAR model.
#' @return The result is returned invisibly if plot is TRUE.
#' @author Matteo Farne', Angela Montanari, \email{matteo.farne2@@unibo.it}
#' @seealso \code{\link[vars]{VAR}}.
#' @examples
#' 	RealGdp.rate.ts<-euro_area_indicators[,1]
#' 	m3.rate.ts<-euro_area_indicators[,2]
#' 	hicp.rate.ts<-euro_area_indicators[,4]	
#' 	cond_bc<-bc_test_cond(RealGdp.rate.ts,m3.rate.ts,hicp.rate.ts,ic.chosen="SC",max.lag=2)
#' @references Breitung, J., Candelon, B., 2006. Testing for short- and long-run causality: A frequency-domain approach. 
#' 	\emph{Journal of Econometrics}. \bold{132}, 2, 363--378.
#' @references Farne', M., Montanari, A., 2018. A bootstrap test to detect prominent Granger-causalities across frequencies. 
#'	<arXiv:1803.00374>, \emph{Submitted}.
#' @export
#' @import vars
#' @importFrom graphics abline par
#' @importFrom stats coef frequency median pf qf quantile residuals spec.pgram
#' @importFrom utils install.packages installed.packages

bc_test_cond<-function(x,y,z,ic.chosen="SC",max.lag=min(4,length(x)-1),plot=F,type.chosen="none",p=0,conf=0.95){

if(length(x)==1){
return("The length of x is only 1")
}

if(length(x)!=length(y)){
return("x and y do not have the same length")
}

if(max.lag>length(x)-1){
return("The chosen number of lags is larger than or equal to the time length")
}

if(!require("vars")){
return("The packages 'vars' could not be found. Please install it to 
proceed.")
}
require(vars)


if(p<=0){
mod=VAR(cbind(x,y,z),lag.max=max.lag,ic=ic.chosen,type.chosen)
if(mod$p>1){
p=mod$p
}
if(mod$p==1){
p=2
mod=VAR(cbind(x,y,z),lag.max=max.lag,ic=ic.chosen,type.chosen)
}
}

if(p>0){
if(p>1){
mod=VAR(cbind(x,y,z),p=p,ic=ic.chosen,type.chosen)
}
if(p==1){
p=2
mod=VAR(cbind(x,y,z),p=p,ic=ic.chosen,type.chosen)
}
}

r<-c(0,0)
r<-matrix(0,2,1)

beta<-matrix(mod$p,1)
for(k in 1:mod$p){
beta[k]<-coef(mod)$x[2+(k-1)*mod$K,1]
}

freq.good=spec.pgram(y,plot=F)$freq/frequency(x);

R<-matrix(0,2,mod$p)
R_all<-array(0,dim=c(2,mod$p,length(freq.good)))
for (l in 1:length(freq.good)){
for(k in 1:mod$p){
R[1,k]<-cos(freq.good[l]*pi*k)
R[2,k]<-sin(freq.good[l]*pi*k)
}
R_all[,,l]<-R;
}

X_design<-cbind(x,y)
n<-dim(X_design)[1]
X_all<-matrix(0,n-as.numeric(mod$p),mod$p)
X_past<-matrix(0,n-as.numeric(mod$p),1)
for(k in 1:mod$p){
X_past<-as.matrix(X_design[seq(1+k,1+k+n-as.numeric(mod$p)-1)-1,2])
X_all[1:(n-as.numeric(mod$p)),k]<-X_past
}
X_all<-X_all[,dim(X_all)[2]:1]

res_yes<-residuals(mod)[,1]
sse<-sum(t(res_yes)%*%res_yes)

F_test<-vector(mode="numeric",length(freq.good))
F_test_pre<-vector(mode="numeric",length(freq.good))
F_thr<-vector(mode="numeric",length(freq.good))
for (l in 1:length(freq.good)){
if(l<length(freq.good)){
r_matrix<-r;
R_matrix<-as.matrix(R_all[,,l]);
F_test_pre[l]<-t(r_matrix-R_matrix%*%beta)%*%solve(R_matrix%*%solve(t(X_all)%*%(X_all))%*%t(R_matrix))%*%(r_matrix-R_matrix%*%beta)
F_test[l]<-(F_test_pre[l]/2)/(sse/(n-2*mod$p))}
if(l==length(freq.good)){
r_matrix<-r[1];
R_matrix<-matrix(R_all[1,,l],1,p);
F_test_pre[l]<-t(r_matrix-R_matrix%*%beta)%*%solve(R_matrix%*%solve(t(X_all)%*%(X_all))%*%t(R_matrix))%*%(r_matrix-R_matrix%*%beta)
F_test[l]<-(F_test_pre[l]/1)/(sse/(n-1*mod$p))
}
}

quant_bc<-vector(mode="numeric",length(freq.good))
for (l in 1:length(freq.good)){
quant_bc[l]<-pf(F_test[l],2,n-2*mod$p,lower.tail=F)##<alpha means significant!
if(l==length(freq.good)){
quant_bc[l]<-pf(F_test[l],1,n-1*mod$p,lower.tail=F)##<alpha means significant!
}
}

alpha<-1-conf
signif_bc<-which(quant_bc<alpha)
length(signif_bc)
F_thr[1:length(freq.good)-1]<-qf(conf,2,n-2*mod$p)
F_thr[length(freq.good)]<-qf(conf,1,n-1*mod$p)

GG<-list(
freq.good,length(x),
freq.good[signif_bc],
conf,
F_test,
F_thr,
roots(mod),
as.numeric(mod$p)
)

names(GG)<-c("frequency","n","confidence_level","significant_frequencies","F-test","F-threshold","roots","delays")

if(plot==F){
return(GG)}

if(plot==T){
plot(F_test,type="l")
abline(F_thr,type="l")
}

}
