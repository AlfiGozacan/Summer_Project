library(rjags)

library(statip)

library(latex2exp)

library(MASS)


# 1. Create model


hierarchical_model <- "

  model {
  
    for( i in 1:N ) {
    
      p[i] <- exp( beta[1] + beta[2] * t[i] ) / ( 1 + exp( beta[1] + beta[2] * t[i] ) )
    
      y[i] ~ dbern( p[i] )
    
    }
    
    beta ~ dmnorm( b_0, invB_0 )
  
  }

"

b_0 <- c(0, 0)

B_0 <- diag( 100, 2, 2 )

invB_0 <- solve( B_0 )

data.bayes <- list( y = y_i
                    , N = N
                    , t = t_i
                    , b_0 = b_0
                    , invB_0 = invB_0 )

model <- jags.model( file = textConnection( hierarchical_model )
                     , data = data.bayes )

adapt( object = model
       , n.iter = 1000 )


# 2. Produce output


big_T = 10^4

n.thin = 1

output <- jags.samples( model = model
                        , variable.names = c("beta")
                        , n.iter = big_T * n.thin
                        , thin = n.thin )


# 3. Plot appropriate sample


names(output)

dim(output$beta)

beta.sample <- output$beta

par(mfrow=c(3, 2))

plot( beta.sample[1,,]
      , type = 'l'
      , main = TeX("Trace Plot of $\\beta_1$")
      , xlab = "Iteration"
      , ylab = TeX("$\\beta_1$")
      , ylim = c(min(beta.sample[1,,]), max(beta.sample[1,,])) )

plot( beta.sample[2,,]
      , type = 'l'
      , main = TeX("Trace Plot of $\\beta_2$")
      , xlab = "Iteration"
      , ylab = TeX("$\\beta_2$")
      , ylim = c(min(beta.sample[2,,]), max(beta.sample[2,,])) )


# 4. Make more diagnostic plots


hist( beta.sample[1,,], prob=TRUE, main=TeX("Density of $\\beta_1$"), xlab=TeX("$\\beta_1$") )

lines( density( beta.sample[1,,] ) )

hist( beta.sample[2,,], prob=TRUE, main=TeX("Density of $\\beta_2$"), xlab=TeX("$\\beta_2$") )

lines( density( beta.sample[2,,] ) )

CDF_detail = 1000

x_lst <- seq( min(beta.sample[1,,]), max(beta.sample[1,,]), length.out=CDF_detail )

y_lst <- rep(NA, CDF_detail)

for(i in 1:CDF_detail) {
  
  y_lst[i] <- mean( beta.sample[1,,] <= x_lst[i] )
  
}

plot( x_lst, y_lst, type='l', main=TeX("CDF of $\\beta_1$"), xlab=TeX("$\\beta_1$"), ylab="CDF" )

x_lst <- seq( min(beta.sample[2,,]), max(beta.sample[2,,]), length.out=CDF_detail )

y_lst <- rep(NA, CDF_detail)

for(i in 1:CDF_detail) {
  
  y_lst[i] <- mean( beta.sample[2,,] <= x_lst[i] )
  
}

plot( x_lst, y_lst, type='l', main=TeX("CDF of $\\beta_2$"), xlab=TeX("$\\beta_2$"), ylab="CDF" )

par(mfrow=c(1, 1))

image(kde2d(beta.sample[1,,], beta.sample[2,,], n=100), col=hcl.colors(100, "YlGn", rev=TRUE), xlab=TeX("$\\beta_1$"), ylab=TeX("$\\beta_2$"))
