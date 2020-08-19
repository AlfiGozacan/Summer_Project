### SGLD + MTM algorithm applied to the logistic regression model

### Created on 18/08/20 18:26 by A.Gozacan

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

u_j <- function(j) {
  
  return( sqrt(j) )
  
}

start_pos <- c(0, 1)

accept.probs <- rep(NA, big_T-1) # A list to contain the acceptance probabilities

minibatch_size = 5

interval_size = 20

epsilon = 1e-3

L = c(20, 20) # Take L[1] sub-samples of size L[2] per iteration



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



# 2.5 Make a posterior distribution function, the one we wish to sample from



posterior <- function( beta, b_0, B_0, batch, predictors ) {
  
  term_1 <- prod( p( predictors, beta ) ^ batch * ( 1 - p( predictors, batch ) ) ^ ( 1 - batch ) )
  
  term_2 <- ( 1 / ( 2*pi * sqrt( det(B_0) ) ) ) * exp( -0.5 * t(beta - b_0) %*% solve(B_0) %*% (beta - b_0) )
  
  return( term_1 * term_2 )
  
}



# 3. Create parameter update function



delta <- function(beta, b_0, B_0, batch, predictors, N, n, epsilon_t, eta_t, t) {
  
  term_1 <- deriv_log_prior( beta, b_0, B_0 )
  
  term_2 <- (N/n) * ( sum_deriv_log_likelihood( batch, predictors, beta ) )
  
  return( ( epsilon_t[t] / 2 ) * ( term_1 + term_2 ) + eta_t[t] )
  
}



# 4. Create the step sizes



#epsilon_t <- a * (b + 1:big_T)^(-gamma) # Step size, decreasing with t

epsilon_t <- rep(epsilon, big_T)



# 5. Create the noise



createNoise <- function(big_T, epsilon_t) {
  
  eta_t <- matrix(NA, 2, big_T)
  
  for(t in 1:big_T) {
    
    eta_t[,t] <- mvrnorm(1, c(0, 0), epsilon_t[t] * id)
    
  }
  
  return(eta_t)
  
}



# 6. Create the chain of beta parameters



runSim <- function( m, n, big_T, t_i, y_i, epsilon_t, u_j, start_pos, N ) { # m = Interval size for "leapfrog" method, n = batch size
  
  beta_t <- matrix(NA, 2, big_T)
  
  beta_t[,1] = start_pos # Starting position
  
  beta_i = matrix(NA, 2, m+1) # For the reps in the intervals
  
  eta_t <- createNoise(big_T, epsilon_t)
  
  var_log_fracs <- rep(NA, big_T-1) # A list to contain the variances of the log acceptance ratios for each iteration
  
  log_fracs <- rep(NA, L[1])
  
  start_time <- Sys.time()
  
  for(t in 1:(big_T-1)) {
    
    batch.indices <- sample(1:N, n, replace=F)
    
    y_batch <- y_i[batch.indices]
    
    t_batch <- t_i[batch.indices]
    
    beta_i[,1] <- beta_t[,t]
    
    post_probs <- rep(NA, m) # The posterior probabilities of the given beta in the chain
    
    for( i in 1:m ) {
      
      beta_i[,i+1] <- beta_i[,i] + delta( beta_i[,i], b_0, B_0, y_batch, t_batch, N, n, epsilon_t, eta_t, t )
      
      post_probs[i] <- posterior( beta_i[,i+1], b_0, B_0, y_batch, t_batch ) * u_j(i)
      
    }
    
    lst <- rmultinom(1, big_T, post_probs) # A list containing the relative frequencies of the posterior probabilities for use in the MTM step
    
    index <- which( lst == max(lst) )[1]
    
    beta_star <- beta_i[,index+1]
    
    proposal <- beta_star
    
    beta_i[,1] <- beta_star
    
    post_probs_reverse <- rep(NA, m) # The posterior probabilities of the given beta in the chain for the second reverse chain
    
    for( i in 1:m ) {
      
      beta_i[,i+1] <- beta_i[,i] + delta( beta_i[,i], b_0, B_0, y_batch, t_batch, N, n, epsilon_t, eta_t, t )
      
      post_probs_reverse[i] <- posterior( beta_i[,i+1], b_0, B_0, y_batch, t_batch ) * u_j(i)
      
    }
    
    ratio <- sum( post_probs ) / sum( post_probs_reverse )
    
    for(l in 1:L[1]) {
      
      batch.indices <- sample(1:N, L[2], replace=F)
      
      y_batch <- y_i[batch.indices]
      
      t_batch <- t_i[batch.indices]
      
      beta_i[,1] <- beta_t[,t]
      
      post_probs <- rep(NA, m) # The posterior probabilities of the given beta in the chain
      
      for( i in 1:m ) {
        
        beta_i[,i+1] <- beta_i[,i] + delta( beta_i[,i], b_0, B_0, y_batch, t_batch, N, L[2], epsilon_t, eta_t, t )
        
        post_probs[i] <- posterior( beta_i[,i+1], b_0, B_0, y_batch, t_batch ) * u_j(i)
        
      }
      
      lst <- rmultinom(1, big_T, post_probs) # A list containing the relative frequencies of the posterior probabilities for use in the MTM step
      
      index <- which( lst == max(lst) )[1]
      
      beta_star <- beta_i[,index+1]
      
      beta_i[,1] <- beta_star
      
      post_probs_reverse <- rep(NA, m) # The posterior probabilities of the given beta in the chain for the second reverse chain
      
      for( i in 1:m ) {
        
        beta_i[,i+1] <- beta_i[,i] + delta( beta_i[,i], b_0, B_0, y_batch, t_batch, N, L[2], epsilon_t, eta_t, t )
        
        post_probs_reverse[i] <- posterior( beta_i[,i+1], b_0, B_0, y_batch, t_batch ) * u_j(i)
        
      }
      
      log_fracs[l] <- log( sum( post_probs ) / sum( post_probs_reverse ) )
      
    }
    
    var_log_fracs[t] <- var(log_fracs)
    
    u <- runif(1, 0, 1)
    
    alpha <- exp ( min( 0, log( ratio ) ) )
    
    accept.probs[t] <- alpha
    
    if( t %% 100 == 0 ){
      
      print(t)
      
    }
    
    if( u < alpha ) {
      
      beta_t[,t+1] <- proposal
      
    }
    
    else {
      
      beta_t[,t+1] <- beta_t[,t]
      
    }
    
  }
  
  end_time <- Sys.time()
  
  run_time <- end_time - start_time
  
  par(mfrow=c(2, 3))
  
  library(latex2exp)
  
  plot(1:big_T, beta_t[1, 1:big_T], type="l", xlab=paste("Iteration"), ylab=TeX("$\\beta_1$"))
  
  abline(h=true_beta[1], col="red")
  
  plot(1:big_T, beta_t[2, 1:big_T], type="l", xlab=paste("Iteration"), ylab=TeX("$\\beta_2$"))
  
  abline(h=true_beta[2], col="red")
  
  image(kde2d(beta_t[1, 1:big_T], beta_t[2, 1:big_T], n=200), col=hcl.colors(100, "YlGn", rev=TRUE), xlab=TeX("$\\beta_1$"), ylab=TeX("$\\beta_2$"))
  
  mtext(paste("SGLD + MTM with Correction:", round(run_time, 3), "mins. Avg acc. prob. =", round(mean(accept.probs), 3)), side = 3, line = -44, outer = TRUE)
  
  return( var_log_fracs )
  
}



# 7. Run the simulation to produce some plots of the convergence



var_log_fracs <- runSim( interval_size, minibatch_size, big_T, t_i, y_i, epsilon_t, u_j, start_pos, N )

par(mfrow=c(1, 1))

plot(1:(big_T-1), var_log_fracs, type='l', xlab="Iteration", ylab=TeX("$s^2$"))





