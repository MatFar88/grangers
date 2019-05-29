#' Conditional Granger-causality estimation
#'
#' \verb{Granger.conditional} calculates the Granger-causality conditional spectrum of a
#' 	time series \verb{x} (effect variable) on a time series \verb{z} (conditioning variable) respect
#'	to a time series \verb{y} (cause variable). It requires package \href{https://CRAN.R-project.org/package=vars}{vars}.
#'
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
#'	Defaults to \verb{''none''}. Alternatives are \verb{c(''none'',''const'',''trend'')}.
#' @param p1 parameter \verb{p} to be passed to function \link[vars]{VAR}.
#'	It corresponds to the number of lags of the first VAR model. Defaults to 0.
#' @param p2 parameter \verb{p} to be passed to function \link[vars]{VAR}.
#'	It corresponds to the number of lags of the second VAR model. Defaults to 0.
#' @description Conditional Granger-causality spectrum was first defined in Geweke (1984). It
#' 	measures the strength of the causal link from time series \verb{y} to time series \verb{x} once
#' 	removed the mediating effect of \verb{z} in the frequency domain. Differently from function
#' 	\code{\link[grangers]{Granger.unconditional}}, this function provides only the unidirectional
#' 	causality from \verb{y} to \verb{x}. Here we need to estimate two VAR models: the first on \verb{x} and \verb{z}, the
#' 	second on \verb{x}, \verb{y}, \verb{z}, by package \href{https://CRAN.R-project.org/package=vars}{vars}. Parameters specified for function VAR hold for
#' 	both estimations. For computational details we refer to Ding et al. (2006).
#' @return \verb{frequency}: frequencies used by Fast Fourier Transform.
#' @return \verb{n}: time series length.
#' @return \verb{Conditional_causality_y.to.x.on.z}:  computed conditional Granger-causality from \verb{y} to \verb{x} on \verb{z}.
#' @return \verb{roots_1}: the roots of the estimated VAR on \verb{x} and \verb{y}.
#' @return \verb{roots_2}: the roots of the estimated VAR on \verb{x}, \verb{y} and \verb{z}.
#' @return The result is returned invisibly if plot is TRUE.
#' @author Matteo Farne', \email{matteo.farne2@@unibo.it}
#' @seealso \code{\link[vars]{VAR}}.
#' @examples
#' 	RealGdp.rate.ts<-euro_area_indicators[,1]
#' 	m3.rate.ts<-euro_area_indicators[,2]
#' 	hicp.rate.ts<-euro_area_indicators[,4]
#' 	cond_m3.to.gdp.by.hicp<-
#'	Granger.conditional(RealGdp.rate.ts,m3.rate.ts,hicp.rate.ts,"SC",4)
#' @references Geweke J., 1984. Measures of conditional linear dependence
#' 	and feedback between time series. \emph{J. Am. Stat. Assoc}. \bold{79}, 907--915.
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

Granger.conditional<-function(x,y,z,ic.chosen="SC",max.lag=min(4,length(x)-1),plot=F,type.chosen="none",p1=0,p2=0){

if(length(x)==1){
return("The length of x is only 1")
}

if(length(x)!=length(y)){
return("x and y do not have the same length")
}

if(length(x)!=length(z)){
return("x and z do not have the same length")
}

if(max.lag>length(x)-1){
return("The chosen number of lags is larger than or equal to the time length")
}

##Ding et al.2006

dd_1<-cbind(x,z)
dd_2<-cbind(x,y,z)

if(!require("vars")){
message("The packages 'vars' could not be found. Please install it to 
proceed.")
}

require(vars)

if (p1==0){
model1=VAR(dd_1,ic=ic.chosen,lag.max=max.lag,type.chosen)
}

if (p1>0){
model1=VAR(dd_1,p=p1,type.chosen)
}

freq.good=spec.pgram(y,plot=F)$freq/frequency(x);

astr_too1<-vector("numeric",length(freq.good));
astr_too2<-vector("numeric",length(freq.good));
astr_too3<-vector("numeric",length(freq.good));
astr_too<-vector("numeric",length(freq.good));
boh_too<-vector("numeric",length(freq.good));

coef_model1=array(0,dim=c(model1$K,model1$K,model1$p))

add<-array(0,dim=c(model1$K,model1$K,model1$p,length(freq.good)))

ADD=array(0,dim=c(model1$K,model1$K,length(freq.good)));

ADD_x1=array(0,dim=c(model1$K,model1$K,length(freq.good)))

H1_1=array(0,dim=c(model1$K,model1$K,length(freq.good)))

if (p2==0){
model2=VAR(dd_2,ic=ic.chosen,lag.max=max.lag,type.chosen)
}

if (p2>0){
model2=VAR(dd_2,p=p2,type.chosen)
}

coef_model2=array(0,dim=c(model2$K,model2$K,model2$p))

add_2<-array(0,dim=c(model2$K,model2$K,model2$p,length(freq.good)))

ADD_2=array(0,dim=c(model2$K,model2$K,length(freq.good)));

ADD_x2=array(0,dim=c(model2$K,model2$K,length(freq.good)));

H1_2=array(0,dim=c(model2$K,model2$K,length(freq.good)));

G_2<-array(0,dim=c(model2$K,model2$K,length(freq.good)));

Q1<-array(0,dim=c(model2$K,model2$K,length(freq.good)));

for (p in 1:model1$K){
for (k in 1:model1$p){
coef_model1[1,p,k]=coef(model1)$x[p+(k-1)*model1$K,1]}}
for (p in 1:model1$K){
for (k in 1:model1$p){
coef_model1[2,p,k]=coef(model1)$z[p+(k-1)*model1$K,1]}}

for (p in 1:model2$K){
for (k in 1:model2$p){
coef_model2[1,p,k]=coef(model2)$x[p+(k-1)*model2$K,1]}}
for (p in 1:model2$K){
for (k in 1:model2$p){
coef_model2[2,p,k]=coef(model2)$y[p+(k-1)*model2$K,1]}}
for (p in 1:model2$K){
for (k in 1:model2$p){
coef_model2[3,p,k]=coef(model2)$z[p+(k-1)*model2$K,1]}}


for (l in 1:length(freq.good)){

for (k in 1:model1$p)
{add[,,k,l]<-coef_model1[,,k]*exp(-2*pi*k*(1i)*freq.good[l])}
}


for (l in 1:length(freq.good)){
for (k in 1:model1$p){
ADD[,,l]=ADD[,,l]+add[,,k,l]}
}

Sigma_model1<-summary(model1)$cov

P_model1<-matrix(c(1,-Sigma_model1[1,2]/Sigma_model1[1,1],0,1),2,2)

for (l in 1:length(freq.good)){
ADD_x1[,,l]=(P_model1)%*%(-ADD[,,l]+diag(x=1,model1$K,model1$K)*1)
}

for (l in 1:length(freq.good)){
H1_1[,,l]=solve(ADD_x1[,,l])
}

summary(model2)

Sigma_model2<-summary(model2)$cov

P1_model2<-matrix(c(1,0,0,-Sigma_model2[2,1]*solve(Sigma_model2[1,1]),1,0,-Sigma_model2[3,1]*solve(Sigma_model2[1,1]),0,1),3,3,byrow=T)
P2_model2<-matrix(c(1,0,0,0,1,0,0,-(Sigma_model2[3,2]-Sigma_model2[3,1]*solve(Sigma_model2[1,1])*Sigma_model2[1,2])*solve(Sigma_model2[2,2]-Sigma_model2[2,1]*solve(Sigma_model2[1,1])*Sigma_model2[1,2]),1),3,3,byrow=T)

P_model2<-P1_model2%*%P2_model2

for (l in 1:length(freq.good)){
for (k in 1:model2$p)
{add_2[,,k,l]<-coef_model2[,,k]*exp(-2*pi*k*(1i)*freq.good[l])}
}


for (l in 1:length(freq.good)){
for (k in 1:model2$p){
ADD_2[,,l]=ADD_2[,,l]+add_2[,,k,l]}
}


for (l in 1:length(freq.good)){
ADD_x2[,,l]=(P_model2)%*%(-ADD_2[,,l]+diag(x=1,model2$K,model2$K)*1)
}

for (l in 1:length(freq.good)){
H1_2[,,l]=solve(ADD_x2[,,l])}

for (l in 1:length(freq.good)){
G_2[1,1,l]<-H1_1[1,1,l]
G_2[1,3,l]<-H1_1[1,2,l]
G_2[3,1,l]<-H1_1[2,1,l]
G_2[3,3,l]<-H1_1[2,2,l]
G_2[2,2,l]<-1
}


for (l in 1:length(freq.good)){
Q1[,,l]<-solve(G_2[,,l])%*%H1_2[,,l]
}


for (l in 1:length(freq.good)){
astr_too1[l]<-Q1[1,1,l]*Sigma_model2[1,1]*Conj(Q1[1,1,l])
astr_too2[l]<-Q1[1,2,l]*Sigma_model2[2,2]*Conj(Q1[1,2,l])
astr_too3[l]<-Q1[1,3,l]*Sigma_model2[3,3]*Conj(Q1[1,3,l])
astr_too[l]<-astr_too1[l]+astr_too2[l]+astr_too3[l]
boh_too[l]<-log(abs(astr_too[l])/abs(astr_too1[l]))
}

Granger_cond_y.to.x.by.z<-boh_too

GG<-list(
freq.good,length(x),
Granger_cond_y.to.x.by.z,
roots(model1),roots(model2)
)


names(GG)<-c("frequency","n","Conditional_causality_y.to.x.on.z","roots_1","roots_2")

if(plot==F){
return(GG)}

if(plot==T){
par(mfrow = c(1, 1))
plot(freq.good,Granger_cond_y.to.x.by.z,type="l",main="Conditional Granger-causality y to x on z")
}

}

