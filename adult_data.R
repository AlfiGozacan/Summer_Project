adult <- read.csv("adult.data", header = FALSE)

adult

N = 10000

t_i <- adult$V1[1:N]/100

hist(t_i)

c(min(t_i), max(t_i))

y_i <- rep(NA, N)

for(i in 1:N) {
  
  if(adult$V15[i] == " <=50K") {
    
    y_i[i] <- 0
    
  } else {
    
    y_i[i] <- 1
    
  }
  
}

hist(y_i)
