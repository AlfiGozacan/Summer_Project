### SGLD + VR algorithm applied to the logistic regression model

### Created on 06/08/20 14:43 by A.Gozacan

### To-do list:

### Bugs:

## The algorithm:



# 0. Some constants



b_0 = c(0, 1) # The expectation of the prior distribution of beta

id = matrix(c(1, 0, 0, 1), 2, 2)

B_0 = id * 100

N = 10000 # The no. of samples in the data set

true_beta = c(0, 1) # Chosen as the true value of beta

t_i <- runif(N, 0, 1)

big_T = 10000 # No. of sweeps

a = 0.017 # These numbers were chosen after trialling different combinations
b = -0.15
gamma = 0.501

start_pos <- c(0, 1)

minibatch_size = 5

interval_size = 20

epsilon = 1e-3



# 1. Produce a synthetic data set



library(MASS)

p <- function( t, beta ) {
  
  exp( beta[1] + beta[2]*t ) / ( 1 + exp( beta[1] + beta[2]*t ) )
  
}

synth_data <- function( b_0, B_0, N, true_beta, t_i ) {
  
  y_i = rep(NA, N)
  
  p_i <- p(t_i, true_beta)
  
  for(i in 1:N) {
    
    y_i[i] <- rbinom(1, 1, p_i[i])
    
  }
  
  return( y_i )
  
}

y_i <- synth_data( b_0, B_0, N, true_beta, t_i )



# 2. Create functions for the log-likelihood and log-prior



deriv_log_prior <- function( beta, b_0, B_0 ) {
  
  B_0inv <- solve(B_0)
  
  vec <- B_0inv %*% (b_0 - beta)
  
  return( c( vec[1], vec[2] ) )
  
}

sum_deriv_log_likelihood <- function( y_i, t_i, beta ) {
  
  beta_1 <- beta[1]
  
  beta_2 <- beta[2]
  
  n = length(y_i)
  
  sum <- c(0, 0)
  
  for( i in 1:n ) {
    
    t <- t_i[i]
    
    y <- y_i[i]
    
    p <- p( t, beta )
    
    scalar <- exp( beta_1 + beta_2 * t ) / ( 1 + exp( beta_1 + beta_2 * t ) ) ^ 2
    
    dp_dbeta <- c( scalar, t * scalar )
    
    sum = sum + dp_dbeta * ( y - p ) / ( p * ( 1 - p ) )
    
  }
  
  return( sum )
  
}



# 3. Create parameter update function



delta <- function(beta, b_0, B_0, batch, predictors, N, n, epsilon_t, g_tilde, beta_tilde, eta_t, t) {
  
  term_1 <- deriv_log_prior( beta, b_0, B_0 )
  
  term_2 <- (N/n) * ( sum_deriv_log_likelihood( batch, predictors, beta )
                        - sum_deriv_log_likelihood( batch, predictors, beta_tilde ) )
  
  term_3 <- g_tilde
  
  return( ( epsilon_t[t] / 2 ) * ( term_1 + term_2 + term_3 ) + eta_t[t] )
  
}



# 4. Create the step sizes



#epsilon_t <- a * (b + 1:big_T)^(-gamma) # Step size, decreasing with t

epsilon_t <- rep(epsilon, big_T)



# 5. Create the noise



createNoise <- function( big_T, epsilon_t ) {
  
  eta_t <- matrix(NA, 2, big_T)
  
  for(t in 1:big_T) {
    
    eta_t[,t] <- mvrnorm(1, c(0, 0), epsilon_t[t] * id)
    
  }
  
  return(eta_t)

}



# 6. Create the chain of beta parameters



runSim <- function( m, n, big_T, t_i, y_i, epsilon_t, start_pos, N ) { # m = Interval size for "leapfrog" method, n = batch size
  
  beta_t <- matrix(NA, 2, big_T)
  
  beta_t[,1] = start_pos # Starting position
  
  beta_tilde <- beta_t[,1]
  
  g_tilde_0 <- sum_deriv_log_likelihood( y_i, t_i, beta_tilde )
  
  eta_t <- createNoise(big_T, epsilon_t)
  
  start_time <- Sys.time()
  
  for(t in 1:(big_T-1)) {
    
    beta_tilde <- beta_t[,1]
    
    g_tilde <- g_tilde_0
    
    if( t %% m == 0 ) {
      
      beta_tilde <- beta_t[,t]
      
      g_tilde <- sum_deriv_log_likelihood( y_i, t_i, beta_tilde )
      
    }
    
    batch.indices <- sample(1:N, n, replace=F)
    
    y_batch <- y_i[batch.indices]
    
    t_batch <- t_i[batch.indices]
      
    beta_t[,t+1] <- beta_t[,t] + delta( beta_t[,t], b_0, B_0, y_batch, t_batch, N, n, epsilon_t, g_tilde, beta_tilde, eta_t, t )
    
  }
  
  end_time <- Sys.time()
  
  run_time <- end_time - start_time
  
  par(mfrow=c(2, 3))
  
  library(latex2exp)
  
  plot(1:big_T, beta_t[1, 1:big_T], type="l", xlab=paste("Iterations"), ylab=TeX("$\\beta_1$"))
  
  abline(h=true_beta[1], col="red")
  
  plot(1:big_T, beta_t[2, 1:big_T], type="l", xlab=paste("Iterations"), ylab=TeX("$\\beta_2$"))
  
  abline(h=true_beta[2], col="red")
  
  image(kde2d(beta_t[1, 1:big_T], beta_t[2, 1:big_T], n=200, lims=c(-0.5, 0.5, 0.5, 1.5)), col=hcl.colors(100, "YlGn", rev=TRUE), xlab=TeX("$\\beta_1$"), ylab=TeX("$\\beta_2$"))

  mtext(paste("SGLD + VR:", round(run_time, 3), "seconds"), side = 3, line = -44, outer = TRUE)
  
}



# 7. Run the simulation to produce some plots of the convergence



runSim( interval_size, minibatch_size, big_T, t_i, y_i, epsilon_t, start_pos, N )








