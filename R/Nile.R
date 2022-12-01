###############################
# Dynamic quantile model
# 1st-order and 2-nd-order quantile Polynomial model 
# Leo Bastos -- leonardo.bastos AT fiocruz DOT br
# Sat May 17 14:38:57 2014

library('GeneralizedHyperbolic')
library('mvtnorm')

source("R/DQMfunctions.R")

# Preparing the data for a first order polynomial model
Nile.dqm.1st <- list(
  y = Nile,
  FF = matrix(1,nr=1),
  GG = matrix(1,nr=1),
  m0=rep(0,1), 
  C0=1e6*diag(1)
)

# Preparing the data for a second order polynomial model (linear growth model)
Nile.dqm.2nd <- list(
  y = Nile,
  FF = matrix(c(1,0),nr=2),
  GG = matrix(c(1,0,1,1),nr=2),
  m0=rep(0,2), 
  C0=1e6*diag(2)
)


mcmc.sample = seq(1000, 2000, by=5)

dqm.1st.out <- DQMmcmc( Nile.dqm.1st, alpha=c(0.1,0.5,0.9), delta=0.9, M=2000)
dqm.2nd.out <- DQMmcmc( Nile.dqm.2nd, alpha=c(0.1,0.5,0.9), delta=0.9, M=2000)

# DIC (experimental)
DICdqm(dqm.1st.out, mcmc.sample, Nile.dqm.1st)
DICdqm(dqm.2nd.out, mcmc.sample, Nile.dqm.2nd)

# Plots
plot.dqm(dqm.1st.out, mcmc.sample, Nile.dqm.1st, ylab="")
plot.dqm(dqm.2nd.out, mcmc.sample, Nile.dqm.2nd, ylab="")



#save(dqm.1st.out, dqm.2nd.out, Nile.dqm.1st, Nile.dqm.2nd, mcmc.sample, file = "nile2.RData")






# separeted plots

aaa <- dqm.1st.out

# Quantile 10%
q = 1

mean.q10 <- apply(aaa[[q]]$eta[mcmc.sample,],2,mean)
IC.q10 <- apply(aaa[[q]]$eta[mcmc.sample,],2,quantile, probs=c(0.025,0.975))

# Quantile 50%
q = 2

median.q50 <- apply(aaa[[q]]$eta[mcmc.sample,],2,median)
IC.q50 <- apply(aaa[[q]]$eta[mcmc.sample,],2,quantile, probs=c(0.025,0.975))

# Quantile 90%
q = 3

median.q90 <- apply(aaa[[q]]$eta[mcmc.sample,],2,median)
IC.q90 <- apply(aaa[[q]]$eta[mcmc.sample,],2,quantile, probs=c(0.025,0.975))



# Plotting

plot(Nile, ylim=range(IC.q10, IC.q90, Nile), lwd=2)
lines(ts(mean.q10,start=1871, end=1970), t="l", col=2, lwd=2, lty=2)
lines(ts(IC.q10[1,], start=1871, end=1970), t="l", col=2, lty=3,lwd=2)
lines(ts(IC.q10[2,], start=1871, end=1970), t="l", col=2, lty=3,lwd=2)



plot(Nile, ylim=range(IC.q10, IC.q90, Nile), lwd=2)
lines(ts(median.q50,start=1871, end=1970), t="l", col=3, lwd=2, lty=2)
lines(ts(IC.q50[1,], start=1871, end=1970), t="l", col=3, lty=3,lwd=2)
lines(ts(IC.q50[2,], start=1871, end=1970), t="l", col=3, lty=3,lwd=2)



plot(Nile, ylim=range(IC.q10, IC.q90, Nile), lwd=2)
lines(ts(median.q90,start=1871, end=1970), t="l", col=4, lwd=2, lty=2)
lines(ts(IC.q90[1,], start=1871, end=1970), t="l", col=4, lty=3,lwd=2)
lines(ts(IC.q90[2,], start=1871, end=1970), t="l", col=4, lty=3,lwd=2)


plot(Nile, ylim=range(IC.q10, IC.q90, Nile), lwd=2)
for(q in c(1,2,3)){
  lines(ts(apply(aaa[[q]]$eta[mcmc.sample,],2,median),start=1871, end=1970), t="l", col=q+1, lwd=2, lty=2)
}
legend("bottomleft",paste("Q",c(0.1,0.5,0.9)), col=2:4, lty=2, lwd=2)














