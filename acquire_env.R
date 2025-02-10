# functions to acquire environmental variable from the reference
acquire_env <- function(mean_env,mean_sd,sample_size){
  
  env_dataset <- data.frame(matrix(NA, nrow = nrow(mean_env)*sample_size, ncol = ncol(mean_env)))
  names(env_dataset) <- names(mean_env)
  
  for (i in 1:nrow(mean_env))
  {
    for (j in 1:ncol(mean_env))
    {
      mean_value <- as.numeric(mean_env[i,j])
      standard_dev <- as.numeric(mean_sd[i,j])
      # using rnorm() to generate data
      env_dataset[(sample_size*i-(sample_size-1)):(sample_size*i),j] <- rnorm(sample_size, mean = mean_value, sd = standard_dev)
    }
  }
  names(env_dataset) <- colnames(mean_env)
  return(env_dataset)
}
# usage example 
setwd("")
mean_env <- read.csv("env_mean.csv",row.names = 1)
mean_sd <- read.csv("env_sd.csv",row.names = 1)
sample_size = 5 
Wang_2023_Rhizosphere_env <- acquire_env(mean_env = mean_env,mean_sd = mean_sd, sample_size = sample_size)
write.csv(Wang_2023_Rhizosphere_env,"Wang_2023_Rhizosphere_ENV.csv")