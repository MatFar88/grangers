#' Unconditional Granger-causality estimation
#'
#'
#' \verb{Granger.unconditional} calculates the Granger-causality unconditional spectrum of
#' 	a time series \verb{x} (effect variable) respect to a time series \verb{y} (cause variable).
#' 	It requireNamespaces package \href{https://CRAN.R-project.org/package=vars}{vars}.
#'
#' @param x univariate time series.
#' @param y  univariate time series (of the same length of \verb{x}).
#' @param ic.chosen estimation method parameter \verb{ic} to be passed to function \link[vars]{VAR} of
#' 	package \href{https://CRAN.R-project.org/package=vars}{vars}. Defaults to ''SC'' (Schwarz criterion). Alternatives are \verb{c(''AIC'',''HQ'',''SC'',''FPE'')}.
#' @param max.lag maximum number of lags \verb{lag.max} to be passed to function \link[vars]{VAR}.
#' 	Defaults to \verb{min(4, length(x) - 1)}.
#' @param plot logical; if TRUE, it returns the plot of unconditional Granger-causality
#' 	spectra on both directions. Defaults to FALSE.
#' @param type.chosen parameter \verb{type} to be passed to function \link[vars]{VAR}.
#'	Defaults to \verb{''none''}. Alternatives are \verb{c(''none'',''const'',''trend'')}.
#' @param p parameter \verb{p} to be passed to function \link[vars]{VAR}.
#'	Defaults to 0.
#' @description Unconditional Granger-causality spectrum was first defined in Geweke (1982).
#'	It measures the strength of the causal link from time series \verb{y} to time series \verb{x} and
#'	viceversa in the frequency domain. It needs to estimate a VAR model on \verb{x} and \verb{y}
#'	by package \href{https://CRAN.R-project.org/package=vars}{vars}. For computational details we refer to Ding et al. (2006).
#' @return \verb{frequency}: frequencies used by Fast Fourier Transform.
#' @return \verb{n}: time series length.
#' @return \verb{Unconditional_causality_y.to.x}: computed unconditional Granger-causality from \verb{y} to \verb{x}.
#' @return \verb{Unconditional_causality_x.to.y}: computed unconditional Granger-causality from \verb{x} to \verb{y}.
#' @return \verb{roots}: the roots of the estimated VAR on \verb{x} and \verb{y}.
#' @return The result is returned invisibly if plot is TRUE.
#' @author Matteo Farne', Angela Montanari, \email{matteo.farne2@@unibo.it}
#' @seealso \link[vars]{VAR}.
#' @examples
#'	RealGdp.rate.ts<-euro_area_indicators[,1]
#'	m3.rate.ts<-euro_area_indicators[,2]
#'	uncond_m3<-Granger.unconditional(RealGdp.rate.ts,m3.rate.ts,"SC",4)
#' @references Geweke, J., 1982. Measurement of linear dependence and feedback between
#'	multiple time series. \emph{J. Am. Stat. Assoc}. \bold{77}, 304--313.
#' @references Ding, M., Chen, Y., Bressler, S.L., 2006. Granger Causality: Basic Theory and
#' 	Application to Neuroscience, Chap.17. \emph{Handbook of Time Series Analysis
#' 	Recent Theoretical Developments and Applications}.
#' @references Farne', M., Montanari, A., 2018. A bootstrap test to detect prominent Granger-causalities across frequencies. 
#'	<arXiv:1803.00374>, \emph{Submitted}.
#' @export
#' @import vars
#' @importFrom graphics abline par
#' @importFrom stats coef frequency median pf qf quantile residuals spec.pgram
#' @importFrom utils install.packages installed.packages

Granger.unconditional<-function(x,y,ic.chosen="SC",max.lag=min(4,length(x)-1),plot=F,type.chosen="none",p=0){

##Ding et al.2006

if(length(x)==1){
return("The length of x is only 1")
}

if(length(x)!=length(y)){
return("x and y do not have the same length")}

if(max.lag>length(x)-1){
return("The chosen number of lags is larger than or equal to the time length")
}

if(!requireNamespace("vars")){
message("The packages 'vars' could not be found. Please install it to 
proceed.")
}

requireNamespace("vars")

if (p==0){
mod=VAR(cbind(x,y),ic=ic.chosen,lag.max=max.lag,type.chosen)
}

if (p>0){
mod=VAR(cbind(x,y),p=p,type.chosen)
}

coef_mod=array(0,dim=c(mod$K,mod$K,mod$p))

freq.good=spec.pgram(y,plot=F)$freq/frequency(x);

add<-array(0,dim=c(mod$K,mod$K,mod$p,length(freq.good)))

ADD=array(0,dim=c(mod$K,mod$K,length(freq.good)));

H1=array(0,dim=c(mod$K,mod$K,length(freq.good)))

ADD_x=array(0,dim=c(mod$K,mod$K,length(freq.good)));

ADD_2=array(0,dim=c(mod$K,mod$K,length(freq.good)));

H1_x=array(0,dim=c(mod$K,mod$K,length(freq.good)))

H1_2=array(0,dim=c(mod$K,mod$K,length(freq.good)))

astr=vector("numeric",length(freq.good));

astr1=vector("numeric",length(freq.good));

astr3=vector("numeric",length(freq.good));

boh=vector("numeric",length(freq.good));

astr_2=vector("numeric",length(freq.good));

astr_2_1=vector("numeric",length(freq.good));

astr_2_3=vector("numeric",length(freq.good));

boh_2=vector("numeric",length(freq.good));

for (p in 1:mod$K){
for (k in 1:mod$p){
coef_mod[1,p,k]=coef(mod)$x[p+(k-1)*mod$K,1]
}
}

for (p in 1:mod$K){
for (k in 1:mod$p){
coef_mod[2,p,k]=coef(mod)$y[p+(k-1)*mod$K,1]
}
}



for (l in 1:length(freq.good)){
for (k in 1:mod$p)
{add[,,k,l]<--coef_mod[,,k]*exp(-2*pi*k*(1i)*freq.good[l])
}
}

for (l in 1:length(freq.good)){
for (k in 1:mod$p){
ADD[,,l]=ADD[,,l]+add[,,k,l]}
H1[,,l]=solve(ADD[,,l]+diag(x=1,mod$K,mod$K)*exp(-2*pi*0*(1i)*freq.good[l]))
}

##

#y-->x

Sigma<-summary(mod)$cov


P_x<-matrix(c(1,-Sigma[1,2]/Sigma[1,1],0,1),2,2)
P_2<-matrix(c(1,0,-Sigma[2,1]/Sigma[2,2],1),2,2)


for (l in 1:length(freq.good)){
ADD_x[,,l]=(P_x)%*%(ADD[,,l]+diag(x=1,mod$K,mod$K)*exp(-2*pi*0*(1i)*freq.good[l]))
ADD_2[,,l]=(P_2)%*%(ADD[,,l]+diag(x=1,mod$K,mod$K)*exp(-2*pi*0*(1i)*freq.good[l]))
}


for (l in 1:length(freq.good)){
H1_x[,,l]=solve(ADD_x[,,l])#ADD[,,l])
H1_2[,,l]=solve(ADD_2[,,l])
}

for (l in 1:length(freq.good)){
astr1[l]=1/4*H1_x[1,1,l]*Sigma[1,1]*Conj(H1_x[1,1,l])
astr3[l]=1/4*(H1[1,2,l]%*%Sigma[2,2])%*%Conj((H1[1,2,l]))   
astr[l]=astr1[l]+astr3[l]
boh[l]=log(abs(astr[l])/abs(astr1[l]))
astr_2_1[l]=1/4*H1_2[2,2,l]*Sigma[2,2]*Conj(H1_2[2,2,l])
astr_2_3[l]=1/4*(H1[2,1,l]%*%Sigma[1,1])%*%Conj((H1[2,1,l]))  
astr_2[l]=astr_2_1[l]+astr_2_3[l]
boh_2[l]=log(abs(astr_2[l])/abs(astr_2_1[l]))
}

Granger_uncond_y.to.x<-boh
Granger_uncond_x.to.y<-boh_2


GG<-list(
freq.good,length(x),
Granger_uncond_y.to.x,
Granger_uncond_x.to.y,
roots(mod)
)


names(GG)<-c("frequency","n","Unconditional_causality_y.to.x","Unconditional_causality_x.to.y","roots")

if(plot==F){
return(GG)}

if(plot==T){
par(mfrow = c(1, 2))
plot(freq.good,Granger_uncond_y.to.x,main="Granger-causality y to x",type="l") 
plot(freq.good,Granger_uncond_x.to.y,main="Granger-causality x to y",type="l") 
}

}
