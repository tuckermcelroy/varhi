### Script for Var-hidim Analysis

## Note: load package "fractal", "expm"

library(fractal)
library(expm)
library(xtable)

setwd("H:\\VARhidim")

# Generic functions
source("specFact.r")
source("specFactmvar.r")
source("VARMAauto.r")
source("autoVARMA.r")
source("VAR2Wold.r")
source("getGUCD.r")
source("spec.ridge.r")

# Specific functions
source("varhi.getaux.r")
source("varhi.build.r")
source("varhi.refine.r")
source("varhi.fitmodel.r")
source("varhi.checkmodel.r")
source("varhi.foreval.r")
source("varhi.forecomp.r")
source("varhi.lasso.r")
source("varhi.genvar.r")
source("varhi.select.r")
 
 
##########################################
###  load QWI data

# from 1948.Q1 to 2017.Q3
gdpur <- read.table("gdpur_data.dat")
unemp <- gdpur[,2]
scale <- 100
gdp <- scale*log(gdpur[,1])

# QWI data is private (no public category), 19 NAICS
# from 1993.Q1 through 2016.Q2; earnings is 1993.Q2 through 2016.Q1
emp <- as.matrix(read.table("qwiEmp.dat"))
hires <- as.matrix(read.table("qwiHires.dat"))
seps <- as.matrix(read.table("qwiSeps.dat"))
jc <- as.matrix(read.table("qwiJC.dat"))
jd <- as.matrix(read.table("qwiJD.dat"))
earn <- as.matrix(read.table("qwiEarn.dat"))

# NAICS:
# 11  	Agriculture, Forestry, Fishing, and Hunting
# 21  	Mining, Quarrying, and Oil and Gas Extraction
# 22		Utilities
# 23  	Construction
# 31-33	Manufacturing
# 42		Wholesale Trade
# 44-45	Retail Trade
# 48-49	Transportation and Warehousing
# 51		Information
# 52		Finance and Insurance
# 53		Real Estate and Rental and Leasing
# 54		Professional, Scientific, and Technical Services
# 55		Management of Companies and Enterprises
# 56		Administrative and Support and Waste Management and Remediation Services
# 61		Educational Services
# 62		Health Care and Social Assistance
# 71		Arts, Entertainment, and Recreation
# 72		Accomodation and Food Services
# 81		Other Services (except Public Administration)
# 92		Public Administration : NOTE this is omitted

## study my.naics codes as auxiliaries, all data in annual growth rates, 1994.Q2 through 2013.Q4
my.names <- c("Agriculture","Mining","Utilities","Construction","Manufacturing","Wholesale","Retail",
	"Transportation","Information","Finance","Real Estate","Professional","Management","Administrative",
	"Educational","Health Care","Arts","Accomodation","Other")
#my.naics <- c(4,5,7,8,9,10,11,12,13)
my.naics <- seq(1,19)

## qwi.x contains all data of interest; common data span is 1993.Q2 through 2016.Q1
gdp <- gdp[182:273]
ur <- unemp[182:273]
emp <- emp[2:93,]
hires <- hires[2:93,]
seps <- seps[2:93,]
jc <- jc[2:93,]
jd <- jd[2:93,]
earn <- earn[1:92,]
qwi.x <- cbind(emp[,my.naics],hires[,my.naics],seps[,my.naics],jc[,my.naics],jd[,my.naics],earn[,my.naics])
qwi.x <- scale*log(qwi.x)
qwi.x <- cbind(gdp,ur,qwi.x)
qwi.x <- ts(diff(qwi.x,lag=4),start=c(1994,2),frequency=4,names=c("gdp","unemp",paste0("emp.",my.names[my.naics]),
	paste0("hires.",my.names[my.naics]),paste0("seps.",my.names[my.naics]),paste0("jc.",my.names[my.naics]),
	paste0("jd.",my.names[my.naics]),paste0("earn.",my.names[my.naics])))
data <- qwi.x
T.total <- dim(data)[1]

## qwi.y is a restricted sample, for use in doing out-of-sample results
my.samp <- T.total - 29		# fitting to pre-GR period, my.samp = 59
#my.samp <- T.total		# fitting to all time points, my.samp = 88
start.date <- c(1994,2)
end.date <- c(start.date[1]+(my.samp-1) %/% 4, start.date[2]+(my.samp-1) %% 4)
if(end.date[2] > 4) { end.date <- c(end.date[1]+1,end.date[2]-4) }
qwi.y <- window(qwi.x,start.date,end.date)
data.sub <- qwi.y
 
## look at aggregate hires and seps
hires.agg <- ts(diff(scale*log(rowSums(hires)),lag=4),start=c(1994,2),frequency=4)
seps.agg <- ts(diff(scale*log(rowSums(seps)),lag=4),start=c(1994,2),frequency=4)
jc.agg <- ts(diff(scale*log(rowSums(jc)),lag=4),start=c(1994,2),frequency=4)
jd.agg <- ts(diff(scale*log(rowSums(jd)),lag=4),start=c(1994,2),frequency=4)


##########################################################################################
####  Analysis, Results, and Evaluation
####  1. Analysis: run the VAR methods to select best auxiliary variables for the core set
####  2. Results: summarize and display, and check model features
####  3. Evaluation: run forecast comparisons

#### Choose settings for QWI data:

###############################################
## Option A: settings for joint UR and GDP
N.core <- c(2,1)
p.mdl <- 2
 
#######################################
## Option B: settings for GDP
N.core <- 1	# GDP
p.mdl <- 4	# best choice for GDP with my.samp = 59

##########################################
## Option C: settings for UR
N.core <- 2		# UR
p.mdl <- 4


#######################################
########### Analysis ##################
# For any of Options A, B, C, run the analysis

thresh <- .001
N.final <- varhi.select(data.sub,N.core,p.mdl,thresh)


########################################
## Option D: take union of best UR and best GDP series

p.mdl <- 4
thresh <- .001

N.core <- 1	# GDP
N.final.gdp <- varhi.select(data.sub,N.core,p.mdl,thresh)

N.core <- 2	# UR
N.final.ur <- varhi.select(data.sub,N.core,p.mdl,thresh)

N.final <- union(N.final.ur,N.final.gdp)
N.core <- c(2,1)
N.final <- c(N.core,sort(setdiff(N.final,N.core)))
 

########################################
## Option E: use VAR-LASSO

N.core <- c(2,1)
p.mdl <- 4
lam <- 1  # all variables
lam <- 10  # all variables
lam <- 50 # all variables
lam <- 60 # all variables
lam <- 65 # all variables
lam <- 66 # no variables
lam <- 70 # no variables
lam <- 100  # no variables
lambdas <- c(rep(0,length(N.core)),rep(lam,ncol(data.sub)-length(N.core)))
fit.lasso <- varhi.lasso(data.sub,p.mdl,lambdas,ywinit=FALSE)
N.final <- fit.lasso[[1]]




###################################
############## Results ############

# results for thresh = .001, p.mdl = 4
# GDP best: N.final <- c(1, 50)
# UR best: N.final <- c(2, 28, 64)
# GDP-UR union best: order UR first, so that GDP has no instant impact on UR
# N.final <- c(2, 1, 28, 50, 64)
# N.core <- c(2,1)
# GDP-UR joint best: N.final <- c(2,1,13)


## summarize best model
print(colnames(data)[setdiff(N.final,N.core)])

############ Plots ###########################
# plot best series

data.final <- data[,N.final]
plot(data.final,xlab="Year",main="")

# plot core series
#plot(data.final[,c(1,2)],main="",xlab="Year",oma.multi=c(4,0,5,0))
#pdf("corePlots.pdf",width=5,height=4)
par(oma=c(2,0,0,0),mar=c(2,4,2,2)+0.1,mfrow=c(2,1),cex.lab=.8)
plot(data.final[,1],main="",xlab="",ylab="UR",
     yaxt="n",xaxt="n")
axis(1,cex.axis=.5)
axis(2,cex.axis=.5)
plot(data.final[,2],main="",xlab="",ylab="GDP",
     yaxt="n",xaxt="n")
axis(1,cex.axis=.5)
axis(2,cex.axis=.5)
mtext(text="Year",side=1,line=1,outer=TRUE)
dev.off()

# plot auxiliary series
#plot(data.final[,-c(1,2)],main="",xlab="Year",oma.multi=c(4,0,5,0))
#pdf("auxPlots.pdf",width=5,height=6)
par(oma=c(2,0,0,0),mar=c(2,4,2,2)+0.1,mfrow=c(3,1),cex.lab=.8)
plot(data.final[,3],main="",xlab="",ylab="Retail.hires",
     yaxt="n",xaxt="n")
axis(1,cex.axis=.5)
axis(2,cex.axis=.5)
plot(data.final[,4],main="",xlab="",ylab="Finance.seps",
     yaxt="n",xaxt="n")
axis(1,cex.axis=.5)
axis(2,cex.axis=.5)
plot(data.final[,5],main="",xlab="",ylab="Manufacturing.jc",
     yaxt="n",xaxt="n")
axis(1,cex.axis=.5)
axis(2,cex.axis=.5)
mtext(text="Year",side=1,line=1,outer=TRUE)
dev.off()

# plot of aggregate hires and seps
#pdf("hiressepsPlots.pdf",width=5,height=3)
par(oma=c(2,0,0,0),mar=c(2,4,2,2)+0.1,mfrow=c(1,1),cex.lab=.8)
plot(hires.agg,main="",xlab="",ylab="",
     yaxt="n",xaxt="n")
lines(seps.agg,lty=3,lwd=2)
axis(1,cex.axis=.5)
axis(2,cex.axis=.5)
mtext(text="Year",side=1,line=1,outer=TRUE)
dev.off()
 
# plot of aggregate jc and jd
#pdf("jcjdPlots.pdf",width=5,height=3)
par(oma=c(2,0,0,0),mar=c(2,4,2,2)+0.1,mfrow=c(1,1),cex.lab=.8)
plot(jc.agg,main="",xlab="",ylab="",
     yaxt="n",xaxt="n")
lines(jd.agg,lty=3,lwd=2)
axis(1,cex.axis=.5)
axis(2,cex.axis=.5)
mtext(text="Year",side=1,line=1,outer=TRUE)
dev.off()

 

############################################
# fit final model with zero constraints

alpha <- .50
mdl <- varhi.fitmodel(data,p.mdl,N.final,alpha,10^(-12))
phi.arry <- mdl[[2]]
var.pred <- mdl[[3]]
print(mdl[[1]])		# Whittle likelihood
N <- length(N.final)
# check stability
A.mat <- diag(N*(p.mdl+1))
A.mat <- A.mat[seq(1,N*p.mdl),-seq(1,N)]
A.mat[seq(1,N),] <- matrix(phi.arry,c(N,N*p.mdl))
Mod(eigen(A.mat)$values)

# check the model
H <- 48
model.output <- varhi.checkmodel(data,phi.arry,var.pred,N.final,TRUE,FALSE,H)
asymp.se <- model.output[[1]]
resids <- ts(model.output[[2]],start=c(1994,2),frequency=4,names=colnames(data)[N.final])

# plot of residual acf
acf(resids,lag.max=H)

# examine coefficients of final model
for(i in 1:p.mdl) { print(xtable(phi.arry[,,i],digits=3)) }
print(xtable(var.pred,digits=3))

# determine tstats
for(i in 1:p.mdl) { print(xtable(phi.arry[,,i]/asymp.se[,,i],digits=3)) }

# fit final model with no zero constraints
fit.yw <- ar.yw(data.final,aic=FALSE,order.max=p.mdl) 
for(i in 1:p.mdl) { print(xtable(fit.yw$ar[i,,],digits=3)) }
print(xtable(fit.yw$var.pred,digits=3))

# coefficients for driving ur
print(xtable(phi.arry[1,,]/asymp.se[1,,],digits=3))
# coefficients for driving gdp
print(xtable(phi.arry[2,,]/asymp.se[2,,],digits=3))

# compute and display impulse response
maxlag <- 20
m <- length(N.final)
u.mat <- getGUCD(var.pred,m)
u.mat <- u.mat[[1]] %*% diag(sqrt(u.mat[[2]]))
wold.coeff <- VAR2Wold(phi.arry,u.mat,maxlag)

#pdf(file="impulsePlot.pdf",width=10,height=10)
par(mfrow=c(m,m),mar=rep(2,4))
for(j in 1:m)
{
	for(k in 1:m)
	{
		title <- paste(colnames(qwi.x)[N.final[k]]," -> ",colnames(qwi.x)[N.final[j]],sep="")
		plot(ts(wold.coeff[j,k,],start=0,frequency=4),ylab="",xlab="Lag",main=title)
		points(ts(wold.coeff[j,k,],start=0,frequency=4),pty=19)
		abline(h=0,col=grey(.8))
	}
}
dev.off() 

# S-VAR format (the A-Model)
u.inv <- solve(u.mat)
print(xtable(u.inv,digits=3))
for(i in 1:p.mdl) { print(xtable(u.inv %*% phi.arry[,,i],digits=3)) }
print(u.inv %*% var.pred %*% t(u.inv))	# check, should be diag(m)


#####################################
#########  Evaluation ###############

data.fore <- varhi.foreval(data,phi.arry,N.final)

## choose one of the core series to evaluate; 
#	my.index is in 1,2,...,length(N.core)
my.index <- 1 
my.index <- 2 

data.comps <- varhi.forecomp(data,data.fore,p.mdl,my.index,my.samp,N.core,N.final,TRUE,FALSE)
data.comps <- varhi.forecomp(data,data.fore,p.mdl,my.index,my.samp,N.core,N.final,FALSE,TRUE)

# Diebold-Mariano test: positive values indicate mvar is superior to uvar
dm.T <- T.total-my.samp
dm.lags <- ceiling(dm.T/2)
forerr.uni <- c(0,diff(data.comps[[2]]^2))
forerr.mvar <- c(0,diff(data.comps[[1]]^2))
dm.diffs <- forerr.uni - forerr.mvar
dm.diffs <- dm.diffs[-seq(1,my.samp)]
dm.acf <- acf(dm.diffs,lag.max=dm.lags)$acf[-1]*var(dm.diffs)
dm.test <- sum(dm.diffs)/sqrt(dm.T*(var(dm.diffs)+2*sum((1-seq(1,dm.lags)/dm.lags)*dm.acf)))
print(c(sum(dm.diffs),dm.test,1-pt(dm.test,df=length(dm.diffs)-1)))

 
