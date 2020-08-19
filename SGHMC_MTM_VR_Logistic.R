### SGHMC + MTM + VR algorithm applied to the logistic regression model

### Created on 13/08/20 15:39 by A.Gozacan

### To-do list:

### Bugs:

## The algorithm:



# 0. Some constants



b_0 = c(0, 1) # The expectation of the prior distribution of beta

id = matrix(c(1, 0, 0, 1), 2, 2)

B_0 = id * 100

N = 20000 # The no. of samples in the data set

true_beta = c(0, 1) # Chosen as the true value of beta

t_i <- runif(N, 0, 1)

M <- id

big_T = 20000 # No. of sweeps

a = 0.017 # These numbers were chosen after trialling different combinations
b = -0.15
gamma = 0.501

C = id

B_hat = matrix(0, 2, 2)

u_j <- function(j) {
  
  return( sqrt(j) )
  
}

start_pos <- c(0, 1)

minibatch_size = 5

interval_size = 5

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



# 3. Create potential energy function



grad_U <- function(beta, b_0, B_0, batch, predictors, beta_tilde, g_tilde) {
  
  n = length(batch)
  
  term_1 <- - (N/n) * ( sum_deriv_log_likelihood( batch, predictors, beta )
                        - sum_deriv_log_likelihood( batch, predictors, beta_tilde ) )
  
  term_2 <- - deriv_log_prior( beta, b_0, B_0 )
  
  return( term_1 + term_2 - g_tilde )
  
}



# 4. Create the step sizes



#epsilon_t <- a * (b + 1:big_T)^(-gamma) # Step size, decreasing with t

epsilon_t <- rep(epsilon, big_T)



# 5. Create the noise function



createNoise <- function(m, t, C, B_hat, epsilon_t ) {
  
  eta_t = mvrnorm(m, c(0, 0), 2 * ( C - B_hat) * epsilon_t[t] )
  
  return( eta_t )
  
}



# 6. Create the chain of beta parameters



runSim <- function( m, n, M, big_T, t_i, y_i, C, B_hat, epsilon_t, u_j, start_pos, L ) { # m = Interval size for "leapfrog" method, n = batch size
                                                                                         #, L = length of sub-samples for variance calculations
  
  beta_t <- matrix(NA, 2, big_T)
  
  beta_t[,1] = start_pos # Starting position
  
  r_t <- mvrnorm(big_T, c(0, 0), M)
  
  beta_i = matrix(NA, 2, m+1) # For the reps in the intervals
  
  r_i = matrix(NA, 2, m+1)
  
  var_log_fracs <- rep(NA, big_T-1)
  
  log_fracs <- rep(NA, L[1])
  
  accept.probs <- rep(NA, big_T-1) # A list to contain the acceptance probabilities
  
  beta_tilde <- beta_t[,1]
  
  g_tilde_0 <- sum_deriv_log_likelihood( y_i, t_i, beta_tilde )
  
  start_time <- Sys.time()
  
  for(t in 1:(big_T-1)) {
    
    beta_tilde <- beta_t[,1]
    
    g_tilde <- g_tilde_0
    
    if(t %% m == 0) {
      
      beta_tilde <- beta_t[,t]
      
      g_tilde <- sum_deriv_log_likelihood( y_i, t_i, beta_tilde )
      
    }
    
    batch.indices <- sample(1:N, n, replace=F)
    
    y_batch <- y_i[batch.indices]
    
    t_batch <- t_i[batch.indices]
    
    beta_i[,1] <- beta_t[,t]
    
    r_i[,1] <- r_t[t,]
    
    noise <- createNoise( m, t, C, B_hat, epsilon_t )
    
    post_probs <- rep(NA, m) # The posterior probabilities of the given beta in the chain
    
    for( i in 1:m ) {
      
      beta_i[,i+1] <- beta_i[,i] + epsilon_t[t] * r_i[,i]
      
      grad_U_tilde <- grad_U( beta_i[,i+1], b_0, B_0, y_batch, t_batch, beta_tilde, g_tilde )
      
      r_i[,i+1] <- r_i[,i] - epsilon_t[t] * grad_U_tilde - epsilon_t[t] * ( C %*% solve(M) %*% r_i[,i] ) + noise[i,]
      
      post_probs[i] <- posterior( beta_i[,i+1], b_0, B_0, y_batch, t_batch ) * u_j(i)
      
    }
    
    lst <- rmultinom(1, big_T, post_probs) # A list containing the relative frequencies of the posterior probabilities for use in the MTM step
    
    index <- which( lst == max(lst) )[1]
    
    beta_star <- beta_i[,index+1]
    
    proposal <- beta_star
    
    beta_i[,1] <- beta_star
    
    r_i[,1] <- r_t[t,]
    
    post_probs_reverse <- rep(NA, m) # The posterior probabilities of the given beta in the chain for the second reverse chain
    
    for( i in 1:m ) {
      
      beta_i[,i+1] <- beta_i[,i] + epsilon_t[t] * - r_i[,i]
      
      grad_U_tilde <- grad_U( beta_i[,i+1], b_0, B_0, y_batch, t_batch, beta_tilde, g_tilde )
      
      r_i[,i+1] <- r_i[,i] - epsilon_t[t] * grad_U_tilde - epsilon_t[t] * ( C %*% solve(M) %*% r_i[,i] ) + noise[i,]
      
      post_probs_reverse[i] <- posterior( beta_i[,i+1], b_0, B_0, y_batch, t_batch ) * u_j(i)
      
    }
    
    ratio <- sum( post_probs ) / sum( post_probs_reverse )
    
    for(l in 1:L[1]) {
      
      batch.indices <- sample(1:N, L[2], replace=F)
      
      y_batch <- y_i[batch.indices]
      
      t_batch <- t_i[batch.indices]
      
      beta_i[,1] <- beta_t[,t]
      
      r_i[,1] <- r_t[t,]
      
      noise <- createNoise( m, t, C, B_hat, epsilon_t )
      
      post_probs <- rep(NA, m) # The posterior probabilities of the given beta in the chain
      
      for( i in 1:m ) {
        
        beta_i[,i+1] <- beta_i[,i] + epsilon_t[t] * r_i[,i]
        
        grad_U_tilde <- grad_U( beta_i[,i+1], b_0, B_0, y_batch, t_batch, beta_tilde, g_tilde )
        
        r_i[,i+1] <- r_i[,i] - epsilon_t[t] * grad_U_tilde - epsilon_t[t] * ( C %*% solve(M) %*% r_i[,i] ) + noise[i,]
        
        post_probs[i] <- posterior( beta_i[,i+1], b_0, B_0, y_batch, t_batch ) * u_j(i)
        
      }
      
      lst <- rmultinom(1, big_T, post_probs) # A list containing the relative frequencies of the posterior probabilities for use in the MTM step
      
      index <- which( lst == max(lst) )[1]
      
      beta_star <- beta_i[,index+1]
      
      beta_i[,1] <- beta_star
      
      r_i[,1] <- r_t[t,]
      
      post_probs_reverse <- rep(NA, m) # The posterior probabilities of the given beta in the chain for the second reverse chain
      
      for( i in 1:m ) {
        
        beta_i[,i+1] <- beta_i[,i] + epsilon_t[t] * - r_i[,i]
        
        grad_U_tilde <- grad_U( beta_i[,i+1], b_0, B_0, y_batch, t_batch, beta_tilde, g_tilde )
        
        r_i[,i+1] <- r_i[,i] - epsilon_t[t] * grad_U_tilde - epsilon_t[t] * ( C %*% solve(M) %*% r_i[,i] ) + noise[i,]
        
        post_probs_reverse[i] <- posterior( beta_i[,i+1], b_0, B_0, y_batch, t_batch ) * u_j(i)
        
      }
      
      log_fracs[l] <- log( sum( post_probs ) / sum( post_probs_reverse ) )
      
    }
    
    var_log_fracs[t] <- var(log_fracs)
    
    u <- runif(1, 0, 1)
    
    alpha <- exp( min( 0, log( ratio ) - 0.5 * var_log_fracs[t] ) )
    
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
  
  plot(1:big_T, beta_t[1, 1:big_T], type="l", xlab=paste("Iterations"), ylab=TeX("$\\beta_1$"))
  
  abline(h=true_beta[1], col="red")
  
  plot(1:big_T, beta_t[2, 1:big_T], type="l", xlab=paste("Iterations"), ylab=TeX("$\\beta_2$"))
  
  abline(h=true_beta[2], col="red")
  
  image(kde2d(beta_t[1, 1:big_T], beta_t[2, 1:big_T], n=200), col=hcl.colors(100, "YlGn", rev=TRUE), xlab=TeX("$\\beta_1$"), ylab=TeX("$\\beta_2$"))
  
  mtext(paste("SGHMC + MTM + VR with Correction:", round(run_time, 3), "mins. Avg acc. prob. =", round(mean(accept.probs), 3)), side = 3, line = -44, outer = TRUE)
  
  return( var_log_fracs )
  
}



# 7. Run the simulation to produce some plots of the convergence



var_log_fracs <- runSim( interval_size, minibatch_size, M, big_T, t_i, y_i, C, B_hat, epsilon_t, u_j, start_pos, L )

par(mfrow=c(1, 1))

plot(1:(big_T-1), var_log_fracs, type='l', xlab="Iteration", ylab=TeX("$s^2$"))







