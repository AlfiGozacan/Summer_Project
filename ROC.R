library(pROC)

n.runs = 200

beta <- runSim( minibatch_size, big_T, t_i, y_i, epsilon_t, start_pos, N )

ROC <- function(n.runs, beta, t_i, y_i, N) {
  
  rates <- matrix(NA, 2, n.runs)
  
  cut_off <- seq(0, 1, length.out=n.runs)

  beta_1 <- mean(beta[1,])

  beta_2 <- mean(beta[2,])
  
  predicted_i <- rep(NA, N)
  
  for(i in 1:N) {
    
    predicted_i[i] <- exp( beta_1 + beta_2 * t_i[i] ) / ( 1 + exp( beta_1 + beta_2 * t_i[i] ) )
    
  }
  
  results <- rep(NA, N)
  
  for(t in 1:n.runs) {
    
    for(i in 1:N) {
      
      if(predicted_i[i] < cut_off[t]) {
        
        results[i] <- 0
        
      } else {
        
        results[i] <- 1
        
      }
      
    }
    
    n.positives <- sum(y_i)
    
    n.negatives <- N - n.positives
    
    conditions <- rep(NA, N)
    
    for(i in 1:N) {
      
      if(results[i] < y_i[i]) {
        
        conditions[i] <- "FN"
        
      }
      
      if(results[i] > y_i[i]) {
        
        conditions[i] <- "FP"
        
      }
      
      if(results[i] == y_i[i] & y_i[i] == 1) {
        
        conditions[i] <- "TP"
        
      }
      
      if(results[i] == y_i[i] & y_i[i] == 0) {
        
        conditions[i] <- "TN"
        
      }
      
    }
    
    n.true_positives <- length(which(conditions=="TP"))
    
    n.false_positives <- length(which(conditions=="FP"))
    
    TPR = n.true_positives / n.positives
    
    FPR = n.false_positives / n.negatives
    
    rates[1,t] = TPR
    
    rates[2,t] = FPR
  
  }
  
  roc <- roc(y_i, predicted_i)
  
  AUC <- auc(roc)
  
  plot(rates[2,]
       , rates[1,]
       , type='l'
       , lwd=2
       , col="Green"
       , xlim=c(0, 1)
       , ylim=c(0, 1)
       , xlab="False Positive Rate"
       , ylab="True Positive Rate"
       , main="ROC Curve for Logistic Regression")
  
  abline(a=0, b=1, lty=2)
  
  text(0.8, 0.2, labels=paste("AUC = ", round(AUC, 4)))
  
}

ROC(n.runs, beta, t_i, y_i, N)
