library(rjags)

library(statip)


# 1. Create model


hierarchical_model <- "

  model {
  
    for( i in 1:N ) {
    
      y[i] ~ dbern( p )
    
    }
    
    p ~ dbeta( a_0, b_0 )
  
  }

"

a_0 = 1

b_0 = 1

data.bayes <- list( y = y_i
                    , N = N
                    , a_0 = a_0
                    , b_0 = b_0 )

model <- jags.model( file = textConnection( hierarchical_model )
                     , data = data.bayes )

adapt( object = model
       , n.iter = 1000 )


# 2. Produce output


big_T = 10000

output <- jags.samples( model = model
                       , variable.names = c("p")
                       , n.iter = big_T )


# 3. Plot appropriate sample


names(output)

dim(output$p)

p.sample <- output$p

p.sample <- p.sample[1,,]

par(mfrow=c(1, 1))

plot( p.sample
      , type = 'l'
      , main = "Trace Plot of p"
      , xlab = "Iteration"
      , ylab = "p"
      , ylim = c(min(p.sample), max(p.sample)) )


# 4. Make more diagnostic plots


a_n = a_0 + N * mean( y_i )

b_n = b_0 + N - N * mean( y_i )

hist( p.sample, prob=TRUE, main="Density of p", xlab="p" )

lines( density( p.sample ) )

curve( dbeta( x, a_n, b_n ), col="Red", add=TRUE )

legend( "topright"
        , legend = c("Empirical KDE", "True Posterior PDF")
        , col = c("Black", "Red")
        , lty=c(1, 1)
        , cex=0.8 )

CDF_detail = 1000

x_lst <- seq( 0, 1, length.out=CDF_detail )

y_lst <- rep(NA, CDF_detail)

for(i in 1:CDF_detail) {
  
  y_lst[i] <- mean( p.sample <= x_lst[i] )
  
}

plot( x_lst, y_lst, type='l', main="CDF of p", xlab="p", ylab="CDF" )

lines( x_lst, pbeta(x_lst, a_n, b_n), col="Red" )

legend( "topright"
        , legend = c("MCMC Approximation", "True Posterior CDF")
        , col = c("Black", "Red")
        , lty=c(1, 1)
        , cex=0.8 )
