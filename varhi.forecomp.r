varhi.forecomp <- function(data,data.fore,p.mdl,my.index,my.samp,N.core,N.final,plot.fores=TRUE,plot.errs=FALSE)
{

	########################################################################
	#	varhi.forecomp by Tucker McElroy, August 15, 2017
	#
	#    Purpose: given a VAR model with final core and auxiliary variables,
	#	compare forecast performance to univariate AR models
 	#    Inputs: 
	#	data: T x N matrix of N time series of length T
	#	data.fore: is T x N.final matrix of 1-step ahead forecasts
	#	p.mdl: selected order of VAR(p.mdl) model
	#	my.index: index of particular core variable of interest
	#	my.samp: index of last time, corresponding to sample for fitting
	#	N.core: indices of core series
	#	N.final: indices of core and final auxiliary series
	#	plot.fores: flag for whether forecasts should be plotted
	#	plot.errs: flag for whether forecast error ratios should be plotted
	#    Outputs:  plots, and
	#	new.errs: forecast errors from the final VAR model
	#	uni.errs: forecast errors from the fitted AR model
	#	uni.fores: forecasts from the fitted AR model
	#
	########################################################################		

T.total <- dim(data)[1]
new.errs <- NULL
for(i in 1:(dim(data.fore)[1]))
{
	new.errs <- c(new.errs, sqrt(sum((data[1:i,N.core[my.index]] - 
		data.fore[1:i,my.index])^2)))
}

data.index <- data[1:my.samp,N.core[my.index]]
mu <- mean(data.index)
data.index <- data.index - mu
fit.uni <- ar.yw(data.index)
 
start.eval <- fit.uni$order
fore.lead <- T.total-start.eval
uni.errs <- NULL
if(fore.lead > 0) {
forecasts <- NULL
for(h in 1:fore.lead)
{
	forecast <- mu  + sum(rev(fit.uni$ar) * (data[(start.eval+h-fit.uni$order):(start.eval+h-1),N.core[my.index]] - mu))
	forecasts <- c(forecasts,forecast)
} 
uni.fores <- c(data[1:start.eval,N.core[my.index]],forecasts) 
for(i in 1:length(uni.fores))
{
	uni.errs <- c(uni.errs, sqrt(sum((data[1:i,N.core[my.index]] - uni.fores[1:i])^2)))
}
}  # end if


########## PLOTS ############################
## red is out-of-sample one-step ahead forecasts (fit is based on my.samp),
#  green is in-sample one-step ahead forecasts

dots <- end(ts(seq(1,my.samp),start=start(data),frequency=frequency(data)))
dots <- dots[1]+(dots[2]-1)/frequency(data)

if(plot.fores) {
plot(ts(data[,N.core[my.index]],start=start(data),frequency=frequency(data)),
	xlab="Year",ylab="",lwd=2)
lines(ts(data.fore[,my.index],start=start(data),frequency=frequency(data)),col=grey(.6),lty=2,lwd=2)
lines(ts(data.fore[1:my.samp,my.index],start=start(data),frequency=frequency(data)),col=grey(.6),lty=5,lwd=2)
lines(ts(uni.fores,start=start(data),frequency=frequency(data)),col=grey(.3),lty=3,lwd=2)
abline(v=dots,lty=1,lwd=2,col=grey(.8))
}

if(plot.errs) {
comp.errs <- ts((new.errs/uni.errs)[-seq(1,20)],start=start(lag(data,-20)),frequency=frequency(data))
#plot(comp.errs,xlab="Year",ylab="",lwd=2,ylim=c(0,1.1))  
plot(comp.errs,xlab="Year",ylab="",lwd=2,ylim=c(min(comp.errs),1.1))  
abline(h=1,lty=2,lwd=2)
abline(v=dots,lty=1,lwd=2,col=grey(.8))
}

return(list(new.errs,uni.errs,uni.fores))
}



