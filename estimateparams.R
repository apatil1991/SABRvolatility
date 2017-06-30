SABRVol<-function(a,b,r,v,F,K,T){
  #--------------------
  #Inputs
  # a = alpha 
  # b = beta
  # r = rho
  # v = vol of vol
  # F = Spot Price
  # K = Strike Price
  # T = Maturity
  #--------------------
  x <- function(z,r){
    log((sqrt(1 - 2*r*z + z^2) + z - r)/(1 - r))
  }
  z <- v/a * (F*K)^((1 - b)/2)*log(F/K)
  Denom <- (F*K)^((1 - b)/2)*
    (1 +
      (1 - b)^2/24*log(F/K)^2 + # 0 if f==K
      (1 - b)^4/1920 *log(F/K)^4 # 0 if f==K
    )
  Term3 <- (1 + ((1 - b)^2/24*a^2/((F*K)^(1 - b)) + (1/4)*(r*b*v*a)/((F*K)^((1 - b)/2)) + (2 - 3*r^2)/24*v^2)*T)
  # To work with vectors of strikes (K's):
  Term2 <- rowSums(cbind(is.na(z/x(z,r)), (1 - is.na(z/x(z,r)))*(z/x(z,r))), na.rm=TRUE)
  # Putting the pieces together and returning the result
  return( (a/Denom) * Term2 * Term3 )
}

Estimateallparams<-function(params, MktStrike, MktVol, F, T, b)
{
  # -------------------------------------------------------------------------
  # Required inputs:
  # MktStrike = Vector of Strikes
  # MktVol = Vector of corresponding volatilities
  # ATMVol = ATM volatility
  # F = spot price
  # T = maturity
  # b = beta parameter
  # -------------------------------------------------------------------------
  a <- params[1]
  r <- params[2]
  v <- params[3]
  N <- length(MktVol)

  ModelVol <- numeric(N)
  error <- numeric(N)
# Define the model volatility and the squared error terms
  for (i in 1:N) {
    ModelVol[i] = SABRVol(a, b, r, v, F, MktStrike[i], T)
    error[i] = (ModelVol[i] - MktVol[i])^2
  }
  
  # Return the SSE
  y <- sum(error, na.rm=TRUE)
  
  # Impose the constraint that -1 <= rho <= +1
  # via a penalty on the objective function
  if (abs(r) > 1 | v<0){
    y <- 1e100
  }
  return(y)
}
