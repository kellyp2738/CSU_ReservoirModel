############################################################################
##
## MODIFICATIONS OF STEVE'S MODEL:
##
## Model of Ebola virus disease (EVD) outbreak in Liberia to compare projections ## made with and without accounting for asymptomatic, non-infectious infections.
##
## Accompanying:
##
## Bellan SE, JRC Pulliam, J Dushoff, and LA Meyers. (2014) Impact of
## asymptomatic infection and acquired immunity on the spread and control of
## Ebola. Commentary, _The Lancet_.
##
## Steve Bellan, September 28, 2014

###############
### LICENSE
###
### This code is made available under a Creative Commons Attribution 4.0
### International License. You are free to reuse this code provided that you
### give appropriate credit, provide a link to the license, and indicate if
### changes were made.
### You may do so in any reasonable manner, but not in any way that suggests
### the licensor endorses you or your use. Giving appropriate credit includes
### citation of the above publication *and* providing a link to this repository:
###
### https://github.com/ICI3D/Ebola
###############

rm(list=ls(all=T))
library(deSolve) ## differential equation solver library
library(Rcpp)
library(RcppEigen)

## CONTROL PARAMETERS
SAVEPLOTS <- FALSE
if(SAVEPLOTS){
	try(system('mkdir Figures'))
	if(file.exists("Figures")){
		plotPath <- 'Figures'
	}else{
		warning("Figures directory does not exist; plotPath will not be defined, and figures will be plotted but not saved.")
	}
}

## Susceptible-Exposed-Infectious/Symptomatic-NonInfectious/Asymptomatic-Recovered Model
## STATE VARIABLES AND DERIVED QUANTITIES:
###   S                  ## susceptible
###   E                  ## exposed (incubating/latent)
###   I                  ## infectious
###   D                  ## dead
###   R                  ## recovered/immune
###   G                  ## gone
###   cumExp             ## cumulative number of exposed/infected individuals (includes asymptomatic)
###   cumInc             ## cumulative number of symptomatic infected individuals (excludes asymptomatic)
###   cumMort            ## cumulative number of EVD deaths
## Set up model equations
seiarModel <- function(t,y,params){
	with(c(as.list(y),params), {
		N <- S + E + I + R ## Total population size
		nu<-mu*N                           ## Births (zero in current model, as are deaths)
		dS <- lambda*R + nu - beta*S*(I+D) - mu*S        ## Susceptible
		dE <- beta*S*(I+D) - sigma*E - mu*E    ## Exposed (incubating)
		dI <- sigma*E - gamma*I - mu*I         ## Infectious
		dD <- cfr*gamma*I - kappa*D    ## Dead but still infectious
		dG <- kappa*D                   ## Gone (dead and no longer infectious)
		dR <- (1-cfr)*gamma*I - mu*R - lambda*R          ## Recovered & immune
		return(list(c(dS=dS,dE=dE,dI=dI,dD=dD,dG=dG,dR=dR)))
	})
}

## The estimated basic reproduction numbers (R0) are
## 1.71 (95% CI, 1.44 to 2.01) for Guinea,
## 1.83 (95% CI, 1.72 to 1.94) for Liberia,
## 2.02 (95% CI, 1.79 to 2.26) for Sierra Leone. (WHO Global Ebola Response Team Lancet 2014)
R0 <- 1.83
N0 <- 4*10^6 ## Liberia's population
N0K <- N0/10^3 ## population size in thousands
## Initialize with one symptomatic infectious individual who will die.
init <- c(S = N0-1,
          E = 0,
          I = 1,
          D = 0,
          G = 0,
          R = 0)
times<-seq(0,500,1) ## Simulate for 500 days

param.vals<-c( ## Other parameters
              beta= 0.02, ## BETA PARAMETER IS CURRENTLY MISCALCULATED Calculated based on R0 and other parameters, see below.
              N0=N0,
              mu= 0.002, ## I MADE THIS NUMBER UP. Assume no birth/death for now, though it doesn't affect this toy model much. For 50 yr life expect mu=.02/365.25
              sigma=1/9.1, ## progression rate = 1/incubation or 1/latent period (assumed to be the
                           ## same for Ebola). Lancet estimat 9.1 days; CDC estimate 6 days.
              lambda=1/4, ## 1/duration immune
              kappa=1/4, ## I MADE THIS NUMBER UP. how long are dead carcasses infectious?
              gamma=1/6,  ## 1/infectious period. CDC estimate 6 days
              cfr = .3) ## case fatality rate. Lancet for Liberia = 72.3%

## Jacobian matrix -- partial derivatives of the SEIDR model
## make.jacobian() fills in paramter values to generate the matrix
make.jacobian<-function(params){
  with(as.list(params),
       matrix(data=c(c(0, 0, -1*beta*init[[1]], -1*beta*init[[1]], lambda),
                c(0, -1*sigma, beta*init[[1]], beta*init[[1]], 0),
                c(0, sigma, gamma, 0, 0),
                c(0, 0, cfr*gamma, kappa, 0),
                c(0, 0, gamma-(cfr*gamma), 0, lambda)), 5, 5)
       )
}

jacobian=make.jacobian(param.vals)
Ro<-max(eigen(jacobian)$values)

## Problems:
## 1. We'd like to set Ro and work backward; this technique doesn't allow that
## 2. The Ro we get here is crazy high -- 93!
## 3. The Ro is not constant with population size (model is density dependent as contacts scale with N)
## 4. Solve the dependence on population size by making the model frequency dependent? (reconsider the partial derivatives)

Ro.calc<-function(params) { ## Analytical solution for R0
    with(as.list(params),
    		 beta*(sigma/(mu+sigma)) * symp * (cfr/(mu+gamma) + (1-cfr)/(mu+gamma))
)}
beta.calc<-function(Ro,params) { ## Solve above function for beta
    with(as.list(params),
    		 Ro/((sigma/(mu+sigma)) * symp * (cfr/(mu+gamma) + (1-cfr)/(mu+gamma)))
)}

## Heffernan et al. 1997 found 14 IgG+ individuals in a population sample of 979 people from 8
## villages in Gabon after an epidemic. We calculate the asymptomatic proportion in the following
## way. The number of symptomatic IgG seropositive post outbreak divided by survival rate (1-CFR)
## gives total # of expected symtomatic infections in this sample (i.e. including those that
## died). The denominator is then symptomatic + asymptomatic infections.
print(
	sympHeffernan <- round((4/.3) / (10 + 4/.3), 2)
)
## Estimated to one significant digit, this gives a rough estimate of 57% for symptomatic proportion in the below
## model example.

## However, some of these 10 asymptomatically infected individuals may have been infected from
## wildlife meaning that the actual asymptomatic infected proportion in this outbreak due to *human*
## transmission would have been smaller. Still, other strong evidence exists for asymptomatic
## infection due to human-to-human transmission. Leroy et al. (2000) found that, of 24 close
## contacts of symptomatic EVD cases, 11 were serpositive but had not shown symptoms. Because we do
## not know how many symptomatic contacts exposed these 24 individuals nor how many close contacts
## these symptomatic individuals had, we cannot calculate the symptomatic proportion. However, we
## can calculate an extremeley conservative lower bound by assuming that

## (1) these 24 individuals corresponded to contacts with 24 symptomatic infected individuals
## (i.e. rather than were multiple contacts for symptomatic infection),
## (2) that the study sampled ALL close contacts of these symptomatic individuals (unlikely)
## (3) R0 = 2.

## Under these conservative assumptions, the 24 symptomatic infections caused 48 symptomatic
## infections, and also 11 asymptomatic infections, yielding a lower bound for the asymptomatic
## proportion of 11/(48+11) = 19%. We suspect the asymptomatic proportion is more likely to be in
## the 50% range based on the above study and due to numerous other studies finding high
## seroprevalence of in non-outbreak populations (see review in Becquart et al. 2010).

####################################################################################################
## Run Model
####################################################################################################

runSEIAR <- function(sympProp, paramVals = param.vals, basicReproNum = R0, browse=F){
	if(browse) browser()
	paramVals['symp'] <- sympProp
	paramVals['beta'] <- beta.calc(basicReproNum,paramVals)
	print(paste("Calculated beta value for ",sympProp*100,"% symptomatic: ",round(paramVals['beta'],3),".",sep=""))
	tc <- data.frame(lsoda(init, times, seiarModel, paramVals))  ## Run ODE model
	tc$N <-  rowSums(tc[,c('S','E','I','D','G','R')])    ## Calculate total population size
	tc[,-1] <- tc[,-1] / 10^3                                    ## Show numbers (other than time) in thousands
	tc$Reff <- R0*(tc$S/tc$N)                                    ## Calculate R_effective
	return(tc)
}

sympVals <- c(1,.5)
tcSymp <- runSEIAR(sympVals[1])     ## 100% symptomatic

## my edit: plotting scripts
plot(tcSymp$time, tcSymp$S, type="l", col="blue")
lines(tcSymp$time, tcSymp$E, col='green')
lines(tcSymp$time, tcSymp$I, col='red')
lines(tcSymp$time, tcSymp$R, col='orange')
lines(tcSymp$time, tcSymp$D, col='black')
lines(tcSymp$time, tcSymp$G, col='yellow')

tail(tcSymp)

tcAsymp <- runSEIAR(sympVals[2])    ## 60% symptomatic
## Compare calculated beta values. Note beta is bigger to make up for lower symptomatic proportion.

tail(tcSymp$N + tcSymp$cumMort) ## Check (population size + cumulative mortality) is constant
tail(tcAsymp$N + tcAsymp$cumMort) ## Check (population size + cumulative mortality) is constant
tshow <- c(1:5,nrow(tcSymp))
show <- c('cumInc','cumExp','cumMort')
tcAsymp[tshow,show]/tcSymp[tshow,show]


## Set calendar time. Meltzer et al. (2014) MMWR estimates 3915 EVD cases in Liberia, Aug 28, 2014 after
## applying a correction factor for unreported EVD cases. We set the day in our model closest to this value
## to be August 28.
aug28 <- tcSymp$time[which.min(abs(tcSymp$cumInc*10^3 - 3915))]
days <- as.Date('2014-08-28') + (tcSymp$time - aug28)
tcSymp$days <- tcAsymp$days <- days ## add calendar days to both modeled time series


## Vaccination threshold calculation
crit <- function(vaccPropNeeded = 1/2, ## required to reduce Reff < 1 for R0 = 2
                 cumIncSymp, ## cumulative prevalence of SYMPTOMATIC cases when vaccination starts
                 cfr=0.7,  ## case fatality rate
                 symp=0.5) { ## symptomatic proportion
	asymp <- cumIncSymp*(1-symp)/symp ## proportion of initial population that are asymptomatic cases (immune)
	surv <- (1-cfr)*cumIncSymp ## number symptomatic cases that survived (immune)
	pop <- 1-cfr*cumIncSymp ## total population size (subtracting dead cases)
	propImmune <- (surv+asymp)/pop ## proportion of population immune
	return(pmax(vaccPropNeeded-propImmune, 0)/(1-propImmune))
}
cI <- 0.15 ## if 10% of population gets EVD
print(crit(cumIncSymp=cI)) ## Need to vaccinate 36% of individuals who have not had EVD if 50% are asymptomatic to reverse epidemic
print(crit(cumIncSymp=cI, symp=1)) ## Need to vaccinate 36% of individuals who have not had EVD if 100% symptomatic

####################################################################################################
## Figures
####################################################################################################

## Figure in Lancet letter
sel <- days > as.Date('2014-09-01') & days < as.Date('2015-01-10')  ## show Sep 2014 - Feb 2015
if(SAVEPLOTS) png(file.path(plotPath, 'rel cumInc 2 panel.png'), w = 4, 5, units='in', res = 300)
    par('ps'=11, mar = c(4.5,5,.5,1), lwd = 2, mgp = c(3,1,0), mfrow = c(2,1))
    ## Comparing cumulative EVD cases with and without accounting for asymptomatic proportion.
    mains <- c("(A) Cumulative # of Cases", '(B) Vaccination Coverage Needed for Elimination')
    mains <- rep('',2)
    ## Percent difference in projected cumulative incidence of symptomatic EVD
    ## cases between symptomatic an asymptomatic models (Panel A)
    ylab <- 'Cumulative Cases\n(Thousands)'
    plot(days[sel], with(tcSymp[sel,], cumInc), lty = 1, bty = 'n', type = 'n', xlab = '', bty='L',
         ylab = ylab, las = 1, main = '', xaxt='n')
    mtext(mains[1], side = 3, line =1)
    lines(days[sel], with(tcSymp[sel,], cumInc), lty = 1)
    lines(days[sel], with(tcAsymp[sel,], cumInc), lty = 2)
    mth <- seq.Date(as.Date('2014-01-01'), as.Date('2015-02-01'), by='month')
    axis.Date(side = 1, at = mth, format='%b %e', las = 2)
    legend('topleft', paste0(sympVals*100, '% Symptomatic'), lty = c(1,2), bty = 'n', cex = 1, bg = 'white')
    cRange <- seq(0, 0.2, by=0.01)
    ## Total vaccination coverage needed (Panel B)
    plot(cRange*N0/10^3, crit(cumIncSymp=cRange), type="l", bty = 'L', las = 1, ylim = c(0,.6), lty = 2, main = '',
         xlab="Cumulative Cases\n(Thousands)", ylab="Target Vaccination \nCoverage")
    mtext(mains[2], side = 3, line=1, at=.09*N0/10^3)
    lines(cRange*N0/10^3, crit(cumIncSymp=cRange, symp=1), lty = 1)
if(SAVEPLOTS) graphics.off()


## Percent difference in projected cumulative incidence of symptomatic EVD
## cases between symptomatic an asymptomatic models
sel <- days > as.Date('2014-10-01') & days < as.Date('2015-01-10')  ## show Oct 2014 - Feb 2015
if(SAVEPLOTS) png(file.path(plotPath, 'rel cumInc.png'), w = 5, 3, units='in', res = 300)
par('ps'=11, mar = c(4,4.5,2,1), lwd = 2, mgp = c(3,1,0))
plot(days[sel], 100*((tcSymp$cumInc/tcAsymp$cumInc)[sel]-1), lty = 1, bty = 'n', type = 'l', xlab = '',
     ylab = '', las = 1, main = '', axes=F)
ticks <- seq(0,50, b = 10)
axis(2, ticks, paste(ticks, '%'), las = 2)
mth <- seq.Date(as.Date('2014-01-01'), as.Date('2015-02-01'), by='month')
axis.Date(side = 1, at = mth, format='%b %e', las = 2)
mtext('Overestimation of Cumulative Incidence of Symptomatic Cases \ndue to Ignoring Asymptomatic Infection',3, cex = 1)
mtext('% Overestimate',2, line = 3.5)
if(SAVEPLOTS) graphics.off()

########################
## Additional Figures ##
########################

## Same as Panel A but on log scale. Shows that trajectory is sub log-linear earlier for asymptomatic
## model.
if(SAVEPLOTS) png(file.path(plotPath, 'cumInc comparison LOG.png'), w = 5, 3, units='in', res = 300)
par('ps'=11, mar = c(4,5.5,.5,.5), lwd = 2, mfcol=c(1,1))
plot(days[sel], tcSymp$cumInc[sel], lty = 1, bty = 'n', type = 'n', xlab = '',  log='y',
     ylab = ylab, las = 1, main = '', xaxt='n')
lines(days[sel], with(tcSymp[sel,], cumInc), lty = 1)
lines(days[sel], with(tcAsymp[sel,], cumInc), lty = 2)
mth <- seq.Date(as.Date('2014-01-01'), as.Date('2015-02-01'), by='month')
axis.Date(side = 1, at = mth, format='%b %e', las = 2)
legend('topleft', paste0(sympVals*100, '% Symptomatic'), lty = 1:2, bty = 'n', cex = 1, bg = 'white')
if(SAVEPLOTS) graphics.off()

## Showing cumulative incidence and Reffective over entire epidemic course
if(SAVEPLOTS) png(file.path(plotPath, 'cumInc & Reff comparison.png'), w = 4, 6, units='in', res = 300)
par(mfrow = c(2,1), 'ps'=12, mar = c(4,5.5,2,.5))
sel <- tcSymp$time < 500
## cumulative incidence
with(tcSymp[sel,], plot(days, cumInc, lty = 1, bty = 'n', type = 'l', xlab = '',
                        ylab = "Cumulative # of Cases\n(thousands)", las = 1, xaxt='n',main = ''))
with(tcAsymp[sel,], lines(days, cumInc, lty = 2))
mth <- seq.Date(as.Date('2014-01-01'), as.Date('2015-12-01'), by='month')
axis.Date(side = 1, at = mth, format='%b %e', las = 2)
legend('topleft', paste0(sympVals*100, '% symptomatic'), lty = 1:2, bty = 'n', cex = .7)
abline(v = which.min(abs(tcSymp$cumInc*10^3 - 3000*2.5)), lty = 3, col = 'green')
## Reff
with(tcSymp[sel,], plot(days, Reff, lty = 1, bty = 'n', type = 'l', xlab = '', xaxt='n',
                        ylab =  expression(R['eff']), las = 1, main=''))
with(tcAsymp[sel,], lines(days, Reff, lty = 2))
axis.Date(side = 1, at = mth, format='%b %e', las = 2)
abline(h=1, lty = 2)
if(SAVEPLOTS) graphics.off()

## Time Series of SEIAR groups
if(SAVEPLOTS) png(file.path(plotPath, 'Role of asymptomatic immunity.png'), w = 7.5, 7, units='in', res = 300)
par(lwd=2, mfcol = c(3,2), 'ps'=12, mar = c(4,5,2,.5))
for(ii in 1:2) {
    modtype <- c('Symp','Asymp')[ii]
    tc <- get(paste0('tc',modtype))
    with(tc, {
        plot(days, E, type="l",xlab='',ylab="thousands", bty = 'n', col = 'orange', las = 1, xaxt='n',
             main = paste0(sympVals[ii]*100, '% symptomatic'), ylim = c(0,1000))
        axis.Date(side = 1, at = mth, format='%b %e', las = 2)
        lines(days,(Isurv + Idead), col="red")
        lines(days,A, col="purple")
        lines(days,R, col="blue")
        legend('topleft', c('incubating', 'symtpomatic', 'asymptomatic', 'immune'),
               col = c('orange', 'red','purple','blue'), lwd = 2, bty = 'n')
        ## cumulative cases
        plot(days, cumInc, type="l",xlab='',ylab="thousands", bty = 'n', col = 'red', las = 1, xaxt='n',
             ylim = c(0,6000))
        lines(days,cumMort, col="brown")
        axis.Date(side = 1, at = mth, format='%b %e', las = 2)
        legend('topleft', c('cumulative incidence\n(symptomatic)', 'cumulative mortality'),
               col = c('red','brown'), lwd = 2, bty = 'n')
        ## Reffective
        plot(days, Reff, type="l",xlab='',ylab=expression(R['eff']), bty = 'n', col = 'red', las = 1, xaxt='n',
             ylim = c(0,R0*1.2), main = expression(R['eff']))
        abline(h=1, lty = 2)
        abline(v=days[which.min(abs(Reff-1))])
        axis.Date(side = 1, at = mth, format='%b %e', las = 2)
    })
}
if(SAVEPLOTS) graphics.off()
