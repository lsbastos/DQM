# Dynamic Quantile Regression
# Leo Bastos
# Functions
# Sat May 17 14:38:57 2014

# First order univariate dynamic quantile model 
DQM1stmcmc.univar <- function(y, alpha=0.5, delta=0.9, M=2000, m0=0, C0=1e6, nsig=0.001, ssig=0.001){
  # Q_\alpha(Y_t) = \theta_t
  
  # Que é equivalente a
  # Y_t | U_t, \theta_t, \sigma, \alpha ~ N(\theta_t + A_alpha U_t, \sigma U_t B_\alpha)
  
  # \theta_t = \theta_{t-1} + w_t, w_t ~ N(0, W_t)
  # \theta_0 ~ N(m0, C0)
  # Discount factor
  # W_t = (1-\delta)/\delta C_{t-1}
  
  # U_t ~ Exp(\sigma^{-1})
  
  # sigma \sim InvGam(nsig/2, ssig/2)
  
  N <- length(y)
  
  # Quantidades auxiliares
  Aalpha <- (1-2*alpha) / (alpha*(1-alpha))
  Balpha <- 2 / (alpha*(1-alpha))
  
  
  # Initial values
  theta.mcmc <- Ut.mcmc <- matrix(nr=M, nc=N)
  theta.mcmc[1,] <- rnorm(N) 
  Ut.mcmc[1,] <- 1
  #theta0.mcmc <- sigma.mcmc  <-  W.mcmc <- numeric(M)
  #sigma.mcmc[1] <- W.mcmc[1] <- 1
  theta0.mcmc <- sigma.mcmc  <-  numeric(M)
  sigma.mcmc[1] <- 10
  
  loglike.mcmc <-  numeric(M)
  
  # auxiliar variables for the FFBS
  mt.ffbs <- numeric(N+1)
  Ct.ffbs <- numeric(N+1)
  
  at.ffbs <- numeric(N)
  Rt.ffbs <- numeric(N)
  
  # Discount factor
  delta.ffbs <- delta
  
  mt.ffbs[1] = m0
  Ct.ffbs[1] = C0
  
  loglike.mcmc[1] <- sum(dnorm(y,theta.mcmc[1,] + Aalpha * Ut.mcmc[1,] , sqrt(sigma.mcmc[1]*Balpha * Ut.mcmc[1,]),log=T))
  
  # MCMC
  nsig.star <- nsig + 3*N
  for( k in 2:M){
    # Sampling sigma
    vec.aux.sigma = (y - theta.mcmc[k-1,] - Aalpha*Ut.mcmc[k-1,])
    
    ssig.star <- ssig + drop(crossprod(vec.aux.sigma/(Ut.mcmc[k-1,]*Balpha),vec.aux.sigma)) + 2*sum(Ut.mcmc[k-1,]) 
    sigma.mcmc[k] <- 1/rgamma(1,nsig.star/2, ssig.star/2)
    
    #  # Sampling W
    #  vec.aux.W = (theta.mcmc[k-1,] - c(theta0.mcmc[k-1],theta.mcmc[k-1,-N]))
    #  
    #  nW.star <- nW + N
    #  sW.star <- sW + drop(crossprod(vec.aux.W))
    #  W.mcmc[k] <- 1/rgamma(1, nW.star/2, sW.star/2)
    
    # Sampling Ut and starting the FFBS
    
    psi.gig <- Aalpha^2 / (Balpha * sigma.mcmc[k]) + 2 / sigma.mcmc[k]
    for(t in 1:N){    
      # Sampling latent variables Ut
      chi.t.gig <- (y[t]-theta.mcmc[k-1,t])^2/(Balpha * sigma.mcmc[k])
      Ut.mcmc[k,t] <- rgig(1, chi.t.gig, psi.gig, 1/2)
      
      # Foward filtering
      # a_t = G_t m_{t-1} .. (m_{t-1}) and R_t = G_t C_{t-1} G_t' + W_t .. (C_t-1 + 1/psi2)
      at.ffbs[t] <- mt.ffbs[t]; Rt.ffbs[t] <- Ct.ffbs[t] / delta.ffbs
      
      # f_t = F_t a_{t} .. (m_{t-1}) and Q_t = F_t R_{t} F_t' + V_t .. (R_t + B_alpha U_t / psi1)
      ft <- at.ffbs[t] + Aalpha * Ut.mcmc[k,t]; Qt <- Rt.ffbs[t] + Balpha * Ut.mcmc[k,t] * sigma.mcmc[k]
      
      # m_t = a_t + R_t F_t' Q_t^{-1} (y_t - f_t) and C_t = R_t - R_t F_t' Q_t^{-1} F_t R_t
      mt.ffbs[t+1] <- at.ffbs[t] + Rt.ffbs[t] / Qt * (y[t]-ft)
      Ct.ffbs[t+1] <- Rt.ffbs[t] - Rt.ffbs[t]^2/Qt
    }
    
    # Backward sampling
    theta.mcmc[k,N] <- rnorm(1, mt.ffbs[N+1], sqrt(Ct.ffbs[N+1])  )
    for(t in (N-1):1){
      # h_t = m_t + C_t G_{t+1}'(R_{t+1})^{-1} (theta_{t+1}-a_{t+1})
      ht <- mt.ffbs[t+1] + Ct.ffbs[t+1] / (Rt.ffbs[t+1]) * (theta.mcmc[k,t+1]-at.ffbs[t+1]) 
      
      # H_t = C_t + C_t G_{t+1}' (R_{t+1})^{-1} G_{t+1}' C_t
      Ht <- Ct.ffbs[t+1] + Ct.ffbs[t+1]^2 / Rt.ffbs[t+1]
      
      theta.mcmc[k,t] <- rnorm(1, ht, sqrt(Ht)  )
    }
    ht <- mt.ffbs[1] + Ct.ffbs[1] / Rt.ffbs[1] * (theta.mcmc[k,2]-at.ffbs[1]) 
    Ht <- Ct.ffbs[1] + Ct.ffbs[1]^2 / Rt.ffbs[1]
    theta0.mcmc[k] <- rnorm(1, ht, sqrt(Ht))
    
    # log(p(y|theta,U,sigma))
    loglike.mcmc[k] <- sum(dnorm(y,theta.mcmc[k,] + Aalpha * Ut.mcmc[k,] , sqrt(sigma.mcmc[k]*Balpha * Ut.mcmc[k,]),log=T))
    
    if(k %% 100 == 0) cat(k,"\t", k/M*100, "% done \n")
    
  }
  
  return( list(theta = cbind(theta0.mcmc,theta.mcmc), Ut = Ut.mcmc, sigma = sigma.mcmc, loglike = loglike.mcmc) )
}


# First order dynamic quantile model 
DQM1stmcmc <- function(y, alpha = c(.25, .5, .75), delta=0.9, M=2000, m0=0, C0=1e6, nsig=0.001, ssig=0.001){
  # Q_\alpha(Y_t) = \theta_t
  
  # Que é equivalente a
  # Y_t | U_t, \theta_t, \sigma, \alpha ~ N(\theta_t + A_alpha U_t, \sigma U_t B_\alpha)
  
  # \theta_t = \theta_{t-1} + w_t, w_t ~ N(0, W_t)
  # \theta_0 ~ N(m0, C0)
  # Discount factor
  # W_t = (1-\delta)/\delta C_{t-1}
  
  # U_t ~ Exp(\sigma^{-1})
  
  # sigma \sim InvGam(nsig/2, ssig/2)
  
  N <- length(y)
  
  # Number of quantiles of interest
  Q <- length(alpha)

  gibbs <- vector(length = Q, mode = "list" )
  names(gibbs) <- paste("q", alpha, sep = "")
    
  # Quantidades auxiliares
  Aalpha <- (1-2*alpha) / (alpha*(1-alpha))
  Balpha <- 2 / (alpha*(1-alpha))
  
  
  # Parameters
  theta.mcmc <- Ut.mcmc <- array( dim = c(Q, M, N))
  
  loglike.mcmc <- theta0.mcmc <- sigma.mcmc  <-  array( dim=c(Q,M))

  
  
  # Initial values
  theta.mcmc[,1,] <- rnorm(N) 
  Ut.mcmc[,1,] <- 1
  
  sigma.mcmc[,1] <- 10
  
  
  # auxiliar variables for the FFBS
  mt.ffbs <- array(dim=c(Q,N+1))
  Ct.ffbs <- array(dim=c(Q,N+1))
  
  at.ffbs <- array(dim=c(Q,N))
  Rt.ffbs <- array(dim=c(Q,N))
  
  # Discount factor
  delta.ffbs <- delta
  
  mt.ffbs[,1] = m0
  Ct.ffbs[,1] = C0
  
  for(q in 1:Q)
    loglike.mcmc[q,1] <- sum(dnorm(y,theta.mcmc[q,1,] + Aalpha[q] * Ut.mcmc[q,1,] , sqrt(sigma.mcmc[q,1]*Balpha[q] * Ut.mcmc[q,1,]),log=T))
  
  # MCMC

  # Variables that are invariant on k and q
  nsig.star <- nsig + 3*N
  
  for( k in 2:M){
    # varying the quantile
    for(q in 1:Q){
      # Sampling sigma
      vec.aux.sigma = (y - theta.mcmc[q,k-1,] - Aalpha[q]*Ut.mcmc[q,k-1,])
      
      ssig.star <- ssig + drop(crossprod(vec.aux.sigma/(Ut.mcmc[q,k-1,]*Balpha[q]),vec.aux.sigma)) + 2*sum(Ut.mcmc[q,k-1,]) 
      sigma.mcmc[q,k] <- 1/rgamma(1,nsig.star/2, ssig.star/2)
      
      psi.gig <- Aalpha[q]^2 / (Balpha[q] * sigma.mcmc[q,k]) + 2 / sigma.mcmc[q,k]
      for(t in 1:N){    
        # Sampling latent variables Ut
        chi.t.gig <- (y[t]-theta.mcmc[q,k-1,t])^2/(Balpha[q] * sigma.mcmc[q,k])
        Ut.mcmc[q,k,t] <- rgig(1, chi.t.gig, psi.gig, 1/2)
        
        # Foward filtering
        # a_t = G_t m_{t-1} .. (m_{t-1}) and R_t = G_t C_{t-1} G_t' + W_t .. (C_t-1 + 1/psi2)
        at.ffbs[q,t] <- mt.ffbs[q,t]; Rt.ffbs[q,t] <- Ct.ffbs[q,t] / delta.ffbs
        
        # f_t = F_t a_{t} .. (m_{t-1}) and Q_t = F_t R_{t} F_t' + V_t .. (R_t + B_alpha U_t / psi1)
        ft <- at.ffbs[q,t] + Aalpha[q] * Ut.mcmc[q,k,t]; Qt <- Rt.ffbs[q,t] + Balpha[q] * Ut.mcmc[q,k,t] * sigma.mcmc[q,k]
        
        # m_t = a_t + R_t F_t' Q_t^{-1} (y_t - f_t) and C_t = R_t - R_t F_t' Q_t^{-1} F_t R_t
        mt.ffbs[q,t+1] <- at.ffbs[q,t] + Rt.ffbs[q,t] / Qt * (y[t]-ft)
        Ct.ffbs[q,t+1] <- Rt.ffbs[q,t] - Rt.ffbs[q,t]^2/Qt
      }
      
      # Backward sampling
      theta.mcmc[q,k,N] <- rnorm(1, mt.ffbs[q,N+1], sqrt(Ct.ffbs[q,N+1])  )
      for(t in (N-1):1){
        # h_t = m_t + C_t G_{t+1}'(R_{t+1})^{-1} (theta_{t+1}-a_{t+1})
        ht <- mt.ffbs[q,t+1] + Ct.ffbs[q,t+1] / (Rt.ffbs[q,t+1]) * (theta.mcmc[q,k,t+1]-at.ffbs[q,t+1]) 
        
        # H_t = C_t + C_t G_{t+1}' (R_{t+1})^{-1} G_{t+1}' C_t
        Ht <- Ct.ffbs[q,t+1] + Ct.ffbs[q,t+1]^2 / Rt.ffbs[q,t+1]
        
        theta.mcmc[q,k,t] <- rnorm(1, ht, sqrt(Ht)  )
      }
      ht <- mt.ffbs[q,1] + Ct.ffbs[q,1] / Rt.ffbs[q,1] * (theta.mcmc[q,k,2]-at.ffbs[q,1]) 
      Ht <- Ct.ffbs[q,1] + Ct.ffbs[q,1]^2 / Rt.ffbs[q,1]
      theta0.mcmc[q,k] <- rnorm(1, ht, sqrt(Ht))
      
      # log(p(y|theta,U,sigma))
      loglike.mcmc[q,k] <- sum(dnorm(y,theta.mcmc[q,k,] + Aalpha[q] * Ut.mcmc[q,k,] , sqrt(sigma.mcmc[q,k]*Balpha[q] * Ut.mcmc[q,k,]),log=T))
      
    }
    
    if(k %% 100 == 0) cat(k,"\t", k/M*100, "% done \n")
    
  }
  
  for(q in 1:Q)
    gibbs[[q]] <- list(theta = cbind(theta0.mcmc[q,],theta.mcmc[q,,]), Ut = Ut.mcmc[q,,], sigma = sigma.mcmc[q,], loglike = loglike.mcmc[q,])
  
  return( gibbs )
}




# # Modelo de quantílico dinâmico de primeira ordem 
# DQM1stmcmc.univar <- function(y, alpha=0.5, delta=0.9, M=2000, m0=0, C0=1e6, nsig=0.001, ssig=0.001){
#   # Q_\alpha(Y_t) = \theta_t
#   
#   # Que é equivalente a
#   # Y_t | U_t, \theta_t, \sigma, \alpha ~ N(\theta_t + A_alpha U_t, \sigma U_t B_\alpha)
#   
#   # \theta_t = \theta_{t-1} + w_t, w_t ~ N(0, W_t)
#   # \theta_0 ~ N(m0, C0)
#   # Discount factor
#   # W_t = (1-\delta)/\delta C_{t-1}
#   
#   # U_t ~ Exp(\sigma^{-1})
#   
#   # sigma \sim InvGam(nsig/2, ssig/2)
#   
#   N <- length(y)
#   
#   # Quantidades auxiliares
#   Aalpha <- (1-2*alpha) / (alpha*(1-alpha))
#   Balpha <- 2 / (alpha*(1-alpha))
#   
#   
#   # Initial values
#   theta.mcmc <- Ut.mcmc <- matrix(nr=M, nc=N)
#   theta.mcmc[1,] <- rnorm(N) 
#   Ut.mcmc[1,] <- 1
#   #theta0.mcmc <- sigma.mcmc  <-  W.mcmc <- numeric(M)
#   #sigma.mcmc[1] <- W.mcmc[1] <- 1
#   theta0.mcmc <- sigma.mcmc  <-  numeric(M)
#   sigma.mcmc[1] <- 10
#   
#   loglike.mcmc <-  numeric(M)
#   
#   # auxiliar variables for the FFBS
#   mt.ffbs <- numeric(N+1)
#   Ct.ffbs <- numeric(N+1)
#   
#   at.ffbs <- numeric(N)
#   Rt.ffbs <- numeric(N)
#   
#   # Discount factor
#   delta.ffbs <- delta
#   
#   mt.ffbs[1] = m0
#   Ct.ffbs[1] = C0
#   
#   loglike.mcmc[1] <- sum(dnorm(y,theta.mcmc[1,] + Aalpha * Ut.mcmc[1,] , sqrt(sigma.mcmc[1]*Balpha * Ut.mcmc[1,]),log=T))
#   
#   # MCMC
#   nsig.star <- nsig + 3*N
#   for( k in 2:M){
#     # Sampling sigma
#     vec.aux.sigma = (y - theta.mcmc[k-1,] - Aalpha*Ut.mcmc[k-1,])
#     
#     ssig.star <- ssig + drop(crossprod(vec.aux.sigma/(Ut.mcmc[k-1,]*Balpha),vec.aux.sigma)) + 2*sum(Ut.mcmc[k-1,]) 
#     sigma.mcmc[k] <- 1/rgamma(1,nsig.star/2, ssig.star/2)
#     
#     #  # Sampling W
#     #  vec.aux.W = (theta.mcmc[k-1,] - c(theta0.mcmc[k-1],theta.mcmc[k-1,-N]))
#     #  
#     #  nW.star <- nW + N
#     #  sW.star <- sW + drop(crossprod(vec.aux.W))
#     #  W.mcmc[k] <- 1/rgamma(1, nW.star/2, sW.star/2)
#     
#     # Sampling Ut and starting the FFBS
#     
#     psi.gig <- Aalpha^2 / (Balpha * sigma.mcmc[k]) + 2 / sigma.mcmc[k]
#     for(t in 1:N){    
#       # Sampling latent variables Ut
#       chi.t.gig <- (y[t]-theta.mcmc[k-1,t])^2/(Balpha * sigma.mcmc[k])
#       Ut.mcmc[k,t] <- rgig(1, chi.t.gig, psi.gig, 1/2)
#       
#       # Foward filtering
#       # a_t = G_t m_{t-1} .. (m_{t-1}) and R_t = G_t C_{t-1} G_t' + W_t .. (C_t-1 + 1/psi2)
#       at.ffbs[t] <- mt.ffbs[t]; Rt.ffbs[t] <- Ct.ffbs[t] / delta.ffbs
#       
#       # f_t = F_t a_{t} .. (m_{t-1}) and Q_t = F_t R_{t} F_t' + V_t .. (R_t + B_alpha U_t / psi1)
#       ft <- at.ffbs[t] + Aalpha * Ut.mcmc[k,t]; Qt <- Rt.ffbs[t] + Balpha * Ut.mcmc[k,t] * sigma.mcmc[k]
#       
#       # m_t = a_t + R_t F_t' Q_t^{-1} (y_t - f_t) and C_t = R_t - R_t F_t' Q_t^{-1} F_t R_t
#       mt.ffbs[t+1] <- at.ffbs[t] + Rt.ffbs[t] / Qt * (y[t]-ft)
#       Ct.ffbs[t+1] <- Rt.ffbs[t] - Rt.ffbs[t]^2/Qt
#     }
#     
#     # Backward sampling
#     theta.mcmc[k,N] <- rnorm(1, mt.ffbs[N+1], sqrt(Ct.ffbs[N+1])  )
#     for(t in (N-1):1){
#       # h_t = m_t + C_t G_{t+1}'(R_{t+1})^{-1} (theta_{t+1}-a_{t+1})
#       ht <- mt.ffbs[t+1] + Ct.ffbs[t+1] / (Rt.ffbs[t+1]) * (theta.mcmc[k,t+1]-at.ffbs[t+1]) 
#       
#       # H_t = C_t + C_t G_{t+1}' (R_{t+1})^{-1} G_{t+1}' C_t
#       Ht <- Ct.ffbs[t+1] + Ct.ffbs[t+1]^2 / Rt.ffbs[t+1]
#       
#       theta.mcmc[k,t] <- rnorm(1, ht, sqrt(Ht)  )
#     }
#     ht <- mt.ffbs[1] + Ct.ffbs[1] / Rt.ffbs[1] * (theta.mcmc[k,2]-at.ffbs[1]) 
#     Ht <- Ct.ffbs[1] + Ct.ffbs[1]^2 / Rt.ffbs[1]
#     theta0.mcmc[k] <- rnorm(1, ht, sqrt(Ht))
#     
#     # log(p(y|theta,U,sigma))
#     loglike.mcmc[k] <- sum(dnorm(y,theta.mcmc[k,] + Aalpha * Ut.mcmc[k,] , sqrt(sigma.mcmc[k]*Balpha * Ut.mcmc[k,]),log=T))
#     
#     if(k %% 100 == 0) cat(k,"\t", k/M*100, "% done \n")
#     
#   }
#   
#   return( list(theta = cbind(theta0.mcmc,theta.mcmc), Ut = Ut.mcmc, sigma = sigma.mcmc, loglike = loglike.mcmc) )
# }


# Dynamic quantile model with dscount factor and F and G terms
DQMmcmc <- function(obj, alpha = c(.25, .5, .75), delta=0.9, M=2000, nsig=0.001, ssig=0.001){
  # obj is a list containg the triple: {y, FF, GG, m0, C0}
    
  # Q_\alpha(Y_t) = \theta_t
  
  # Que é equivalente a
  # Y_t | U_t, \theta_t, \sigma, \alpha ~ N(\theta_t + A_alpha U_t, \sigma U_t B_\alpha)
  
  # \theta_t = \theta_{t-1} + w_t, w_t ~ N(0, W_t)
  # \theta_0 ~ N(m0, C0)
  # Discount factor
  # W_t = (1-\delta)/\delta C_{t-1}
  
  # U_t ~ Exp(\sigma^{-1})
  
  # sigma \sim InvGam(nsig/2, ssig/2)
  
  N <- length(obj$y)
  P <- length(obj$FF)
  
  # Number of quantiles of interest
  Q <- length(alpha)
  
  
  
  gibbs <- vector(length = Q, mode = "list" )
  names(gibbs) <- paste("q", alpha, sep = "")
  
  # Quantidades auxiliares
  Aalpha <- (1-2*alpha) / (alpha*(1-alpha))
  Balpha <- 2 / (alpha*(1-alpha))
  
  
  # Parameters
  theta.mcmc <- array( dim = c(Q, M, P, N))
  Ut.mcmc <- array( dim = c(Q, M, N))
  eta.mcmc <- array( dim = c(Q, M, N))
  
  loglike.mcmc <- sigma.mcmc  <-  array( dim=c(Q,M))
  
  theta0.mcmc <-  array( dim=c(Q,M,P))
  
  
  # Initial values
  for(q in 1:Q){
    for(p in 1:P){
      theta.mcmc[q,1,p,] <- rnorm(N,0)
    }
    eta.mcmc[q,1,] <- drop(t(obj$FF) %*% theta.mcmc[q,1,,])
  }
  
  
  
  Ut.mcmc[,1,] <- 1
  sigma.mcmc[,1] <- 10
  
  
  # auxiliar variables for the FFBS
  mt.ffbs <- array(dim=c(Q,P,N+1))
  Ct.ffbs <- array(dim=c(Q,P,P,N+1))
  
  at.ffbs <- array(dim=c(Q,P,N))
  Rt.ffbs <- array(dim=c(Q,P,P,N))
  
  # Discount factor
  delta.ffbs <- delta
  
  
  for(q in 1:Q){
    mt.ffbs[q,,1] = obj$m0
    Ct.ffbs[q,,,1] = obj$C0
    
    loglike.mcmc[q,1] <- sum(dnorm(obj$y, eta.mcmc[q,1,] + Aalpha[q] * Ut.mcmc[q,1,] , sqrt(sigma.mcmc[q,1]*Balpha[q] * Ut.mcmc[q,1,]),log=T))
  }
  
  
  
  # MCMC
  
  # Variables that are invariant on k and q
  nsig.star <- nsig + 3*N
  
  for( k in 2:M){
    # varying the quantile
    for(q in 1:Q){
      # Sampling sigma
      vec.aux.sigma = (obj$y - eta.mcmc[q,k-1,] - Aalpha[q]*Ut.mcmc[q,k-1,])/sqrt(Ut.mcmc[q,k-1,]*Balpha[q])
      
      ssig.star <- ssig + drop(crossprod(vec.aux.sigma)) + 2*sum(Ut.mcmc[q,k-1,]) 
      sigma.mcmc[q,k] <- 1/rgamma(1,nsig.star/2, ssig.star/2)
      
      psi.gig <- Aalpha[q]^2 / (Balpha[q] * sigma.mcmc[q,k]) + 2 / sigma.mcmc[q,k]
      for(t in 1:N){    
        # Sampling latent variables Ut
        chi.t.gig <- (obj$y[t]-eta.mcmc[q,k-1,t])^2/(Balpha[q] * sigma.mcmc[q,k])
        Ut.mcmc[q,k,t] <- rgig(1, chi.t.gig, psi.gig, 1/2)
        
        # Foward filtering
        # a_t = G m_{t-1} and R_t = G C_{t-1} G' + W_t 
        at.ffbs[q,,t] <- obj$GG %*% mt.ffbs[q,,t]; 
        Rt.ffbs[q,,,t] <- tcrossprod(obj$GG %*% Ct.ffbs[q,,,t], obj$GG) / delta.ffbs
        
        # f_t = F'a_{t} and Q_t = F' R_{t} F + V_t 
        ft <- drop(crossprod(obj$FF, at.ffbs[q,,t]) + Aalpha[q] * Ut.mcmc[q,k,t])
        Qt <- drop(crossprod(obj$FF, Rt.ffbs[q,,,t] %*% obj$FF) + Balpha[q] * Ut.mcmc[q,k,t] * sigma.mcmc[q,k])
        
        # m_t = a_t + R_t F' Q_t^{-1} (y_t - f_t) and C_t = R_t - R_t F' Q_t^{-1} F R_t
        At <- Rt.ffbs[q,,,t] %*% obj$FF  / Qt
        mt.ffbs[q,,t+1] <- at.ffbs[q,,t] + At * (obj$y[t]-ft)
        Ct.ffbs[q,,,t+1] <- Rt.ffbs[q,,,t] - tcrossprod(At)  * Qt
      }
      
      # Backward sampling
      theta.mcmc[q,k,,N] <- rmvnorm(1, mean=mt.ffbs[q,,N+1], sigma=as.matrix(Ct.ffbs[q,,,N+1])  )
      for(t in (N-1):1){
        aux.BS <- tcrossprod(Ct.ffbs[q,,,t+1], obj$GG) %*% solve(Rt.ffbs[q,,,t+1])
        
        # h_t = m_t + C_t G_{t+1}'(R_{t+1})^{-1} (theta_{t+1}-a_{t+1})
        ht <- mt.ffbs[q,,t+1] + aux.BS %*% (theta.mcmc[q,k,,t+1]-at.ffbs[q,,t+1]) 
        
        # H_t = C_t + C_t G_{t+1}' (R_{t+1})^{-1} G_{t+1}' C_t
        Ht <- Ct.ffbs[q,,,t+1] - aux.BS %*% obj$GG %*% Ct.ffbs[q,,,t+1]
        
        theta.mcmc[q,k,,t] <- rmvnorm(1, mean=ht, sigma = Ht, method="svd"  )
      }
      
      eta.mcmc[q,k,] <- drop(t(obj$FF) %*% theta.mcmc[q,k,,])
      
      aux.BS <- tcrossprod(Ct.ffbs[q,,,1], obj$GG) %*% solve(Rt.ffbs[q,,,1])
      ht <- mt.ffbs[q,,1] + aux.BS %*% (theta.mcmc[q,k,,1]-at.ffbs[q,,1]) 
      
      # H_t = C_t + C_t G_{t+1}' (R_{t+1})^{-1} G_{t+1}' C_t
      Ht <- Ct.ffbs[q,,,1] - aux.BS %*% obj$GG %*% Ct.ffbs[q,,,1]
      
      theta0.mcmc[q,k,] <- rmvnorm(1, mean=ht, sigma=Ht, method="svd")
      
      # log(p(y|theta,U,sigma))
      loglike.mcmc[q,k] <- sum(dnorm(obj$y, eta.mcmc[q,k,] + Aalpha[q] * Ut.mcmc[q,k,] , sqrt(sigma.mcmc[q,k]*Balpha[q] * Ut.mcmc[q,k,]),log=T))      
    }
    
    if(k %% 100 == 0) cat(k,"\t", k/M*100, "% done \n")
    
  }
  
  # theta.mcmc[1,,1,31]
  # gibbs$q0.25$theta[,1,31]
  
  for(q in 1:Q){
    gibbs[[q]] <- list(theta = theta.mcmc[q,,,], theta0 = theta0.mcmc[q,,], 
                       Ut = Ut.mcmc[q,,], sigma = sigma.mcmc[q,], 
                       loglike = loglike.mcmc[q,], eta = eta.mcmc[q,,]  )
  }
  
  gibbs$quantiles <- alpha
  
  return( gibbs )
}


plot.dqm <- function(output, mcmc.index, obj,...){
  
  Q <- length(output$quantiles)
  
  aux <- vector(Q, mode="list")
  
  range.plot <- range(obj$y)
  
  for(q in 1:Q){
    aux[[q]]$mean <- apply(output[[q]]$eta[mcmc.index,],2,mean)
    aux[[q]]$IC <- apply(output[[q]]$eta[mcmc.index,],2, quantile, probs = c(0.025,0.975))
    range.plot <- range(range.plot, aux[[q]]$IC)
  }
  
  if(is.ts(obj$y)==FALSE)
    obj$y <- as.ts(obj$y)
  
  aux.ts <- tsp(obj$y)
  
  plot(obj$y, ylim=range.plot, lwd=2,...)
  for(q in 1:Q){
    lines(ts(aux[[q]]$mean, start=aux.ts[1], end=aux.ts[2], frequency=aux.ts[3]), col=q+1, lty=2, lwd=2)
    lines(ts(aux[[q]]$IC[1,], start=aux.ts[1], end=aux.ts[2], frequency=aux.ts[3]), col=q+1, lty=3, lwd=2)
    lines(ts(aux[[q]]$IC[2,], start=aux.ts[1], end=aux.ts[2], frequency=aux.ts[3]), col=q+1, lty=3, lwd=2)    
  }
  NULL
}


# Calculating DIC (experimental)
DICdqm <- function(output, mcmc.index, obj){
  
  alpha <- output$quantiles
  
  Aalpha <- (1-2*alpha)/(alpha*(1-alpha))
  
  Balpha <- 2/(alpha*(1-alpha))
  
  Q <- length(alpha)
  
  DIC <- data.frame(DIC=numeric(Q), pD = NA)
  
  rownames(DIC) <- names(output)[-(Q+1)]
  
  for(q in 1:Q){
    # Information criteria
    
    Dbar <- mean(-2*output[[q]]$loglike[mcmc.index])
    
    Dthetabar <- -2*sum(dnorm(obj$y, apply(output[[q]]$eta[mcmc.index,],2,mean) + Aalpha[q] * apply(output[[q]]$Ut[mcmc.index,],2,mean) , 
                              sqrt(mean(output[[q]]$sigma[mcmc.index])*Balpha[q] * apply(output[[q]]$Ut[mcmc.index,],2,"mean")),log=T))
    
    # p_d according to Spigehalther et al. (2002)
    DIC[q,2] <- Dbar - Dthetabar
    
    ## According to Gelman et al. (2004, p. 182)
    # pD <- 0.5 * var(-2*loglike.mcmc[mcmc.sample])
    
    # DIC (Spieghalter et al. 2002)
    DIC[q,1] <- DIC[q,2] + Dbar
    
  }
  
  return(DIC)
}