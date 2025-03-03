rNeymanA <- function(n, lambda, theta){
    r <- numeric()
    for (j in 1:n) {
        k <- rpois(1,lambda)
        r[j] <- sum(rpois(k,theta)) 
    }
    return(r)
}

dNeyman <- function(x, lambda, theta){
 p <- numeric()
 p[1]<- exp(-lambda + lambda*exp(-theta))
 c <- lambda*theta*exp(-theta)
 if(x == 0){
   return(p[1])
 }
 else{
   for (i in 1:x) {
     suma = 0
     for (r in 0:(i-1)) {
       suma = suma + (theta^r)/(factorial(r))*p[i-r]
     }
     p[i+1] = (c/(i))*suma # +1 per l'R
   }
   res <- p[i+1]
   return(res)
 }
}

emNeyman <- function(X, itermax = 500, EPSILON = 1e-5) {
    esAvg <- mean(X)
    esVar <- var(X)
    lambda <- esAvg^2 / (esVar - esAvg)
    theta <- (esVar - esAvg) / esAvg
    converged <- FALSE

    for (iter in 1:itermax) {
        param <- c(esAvg, esVar)

        # E
        temp <- lambda * exp(-theta)
        delta <- temp * (1 + digamma(temp + X) - digamma(temp))

        # M
        lambda <- mean(delta)
        theta <- sum(X) / sum(delta)

        esAvg <- lambda * theta
        esVar <- lambda * theta * (1 + theta)
        param_new <- c(esAvg, esVar)

        # check
        D <- sum((param_new - param)^2)
        if (D < EPSILON) {
            converged <- TRUE
            break
        }
    }

    cat("converged :", converged, "\n")
    cat("lambda :", lambda, "\n")
    cat("theta :", theta, "\n")
}

X <- rNeymanA(1000, 5, 1)
emNeyman(X)

lambda <- 10
theta <- 1
y <- 4

temp <- lambda * exp(- theta)

m <- 1:500
left <- sum(temp^m * m^y / factorial(m))
left
right <- exp(temp) * gamma(temp + y) / gamma(temp) 
right

