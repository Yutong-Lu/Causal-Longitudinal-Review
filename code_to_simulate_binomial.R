# Function to create binomial distributions

# Main function
to_prob <- function(x){
  x <- exp(x)
  x <- x / (x+1)
  return(x)
}

# This converts a value on the log-odds scale to probability

# Log-odds = 0 is probability 0.5
to_prob(0)

# Log-odds = 2 is ~ 0.89
to_prob(2)

# log-odds = -2 is ~ .12
to_prob(-2)

# To add in additional terms on the OR scale, wrap them in log(). OR <1 decreases probability
to_prob(0 + log(0.5))
#Should be prob < 0.5

# OR > 1 increases the probability
to_prob(0 + log(1.5))