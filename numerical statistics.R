###########################################################################################################
#This program implements the maximum likelihood estimation method of a two parameters weibull distribution.
#author: THOMÁS FREUD DE MORAIS GONÇALVES
#date: 09/24/2020
#last revision: 09/24/2020
###########################################################################################################

#in order to run this program we gonna need the following packages.
library(moments)
library(MASS)
library("car")
library("dplyr")

#this function implements a pseudorandom number generator applying the mixed congruential generator model.
randGen = function(size, seed) {
  M = (2 ^ 32) - 1
  a = 1664525
  c = 1013904223
  x0 = seed
  x = c()
  for (i in 1:size) {
    x[i] = (x0 = (a * x0 + c) %% M)
  }
  
  return(x = x / M)
}

#weibull random variables generator using inverse method
rrweibull = function(c, b, size, uniform) {
  memory = c()
  memory = b * (-log(uniform)) ^ (1 / c)
  
}

#log-likelihood vectorized function
f = function(param, x) {
  c = param[1]
  b = param[2]
  n = length(x)
  - sum(log(c / b) + (c - 1) * (log(x / b)) - ((x / b) ^ c))
  
}

#log-likelihood gradient function
df1 = function(param, x) {
  c = param[1]
  b = param[2]
  n = length(x)
  
  c(-sum((1 / c) + log(x) - log(b) - log(x / b) * (x / b) ^ c),-sum((-c / b) + (c *
                                                                                  x ^ c) / b ^ (c + 1)))
  
}
#We gonna save the results in the following vectors
reg_skewness_shape = c()
reg_skewness_scale = c()
reg_kurtosis_shape = c()
reg_kurtosis_scale = c()
reg_sample_size = c()
reg_time_execution = c()
reg_variance_shape = c()
reg_variance_scale = c()
reg_bias_shape = c()
reg_bias_scale = c()
mean_shape = c()
mean_scale = c()

#setting the size of samples for each running of the program
sample_numbers = c(10, 15, 30, 60)

for (hh in 1:length(sample_numbers)) {
  shape = 3
  scale = 3
  dados = c()
  sample_size = sample_numbers[hh]
  
  #getting the start time in order to evaluate how long will take the program to run.
  start_time <- Sys.time()
  
  #number of Monte Carlo replications
  nRep = 10000
  seed = 1990
  
  #Generating a nRep*sample_size size vector of uniform random numbers necessary for generate Weibull samples
  uniform = randGen(sample_size * nRep + 1, seed)
  
  #Generating weibull sample
  x = rrweibull(shape, scale, sample_size * nRep + 1, uniform)
  i = 0
  j = 0
  
  #setting up a progress bar
  pb = txtProgressBar(min = 0, max = nRep, style = 3)
  
  #Here is the optimization program
  while (i < nRep) {
    Sys.sleep(0.0001)
    offset = sample_size * (j)
    j = j + 1
    # par= c(c,b) is the initial guess.
    res = optim(
      par = c(5, 1),
      fn = f,
      gr = df1,
      method = "L-BFGS-B",
      x = x[offset:(1 + offset + sample_size)],
      lower = c(0.01, 0.01)
    )
    
    if (res$convergence != 0) {
      print("Convergence Fail")
      M = (2 ^ 32) - 1
      rm(x) #removing data
      seed = uniform[1 + offset] * M
      uniform = randGen(sample_size * (nRep - i), seed)
      x = rrweibull(shape, scale, sample_size * (nRep - i), uniform)
      i = i - 1
      j = 0
      
    }
    dados = rbind(dados, res$par)
    i = i + 1
    #setting a progress bar
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  end_time <- Sys.time()
  print("TIME OF EXECUTION")
  print(end_time - start_time)
  reg_time_execution[hh] = end_time - start_time
  reg_shape = dados[, 1]
  reg_scale = dados[, 2]
  reg_skewness_shape[hh] = skewness(reg_shape)
  reg_skewness_scale[hh] = skewness(reg_scale)
  reg_kurtosis_shape[hh] = kurtosis(reg_shape)
  reg_kurtosis_scale[hh] = kurtosis(reg_scale)
  reg_sample_size[hh] = sample_size
  reg_variance_shape[hh] = var(reg_shape)
  reg_variance_scale[hh] = var(reg_scale)
  mean_shape[hh] = mean(reg_shape)
  mean_scale[hh] = mean(reg_scale)
  reg_bias_shape[hh] = 100 * (mean_shape[hh] - shape) / shape
  reg_bias_scale[hh] = 100 * (mean_scale[hh] - scale) / scale
  
  
}

#generating tables from the simulated data about time of execution
data_time_execution = data.frame(cbind(reg_sample_size, reg_time_execution))
names(data_time_execution) = c("Samples", "Time of Execution")

#rounding numbers for two decimal places
data_time_execution = data_time_execution %>% mutate_at(vars("Time of Execution"), funs(round(., 2)))

#generating tables from the simulated data
records = data.frame(
  reg_sample_size,
  mean_shape,
  mean_scale,
  reg_skewness_shape,
  reg_skewness_scale,
  reg_kurtosis_shape,
  reg_kurtosis_scale,
  reg_variance_shape,
  reg_variance_scale,
  reg_bias_shape,
  reg_bias_scale
)

names(records) = c(
  "Samples",
  "Shape",
  "Scale",
  "Skewness.Shape",
  "Skewness.Scale",
  "Kurt.Shape",
  "Kurt.Scale",
  "Var.Shape",
  "Var.Scale",
  "Bias.Shape",
  "Bias.Scale"
)

#rounding numbers for two decimal places
records = records %>% mutate_at(vars(-Samples), funs(round(., 2)))

#printing simulated results
print(records)
print("Var = Variance; Kurt = Kurtosis, ")
print(data_time_execution)

#saving data in files
write.csv(data_time_execution, "data_time_execution.txt")
write.csv(records, "records.txt")



###############################################
#From now on we generate some plots of the data.
par(mfrow = c(2, 2))
plot(
  reg_sample_size,
  cumsum(reg_skewness_shape) / (1:length(sample_numbers)),
  type = 'l',
  xlab = "Sample Size",
  ylab = "c - Skewness",
  col = 'blue',
  lwd = 2
)
plot(
  reg_sample_size,
  cumsum(reg_skewness_scale) / (1:length(sample_numbers)),
  type = 'l',
  xlab = "Sample Size",
  ylab = "b - Skewness",
  col = 'red',
  lwd = 2
)
plot(
  reg_sample_size,
  cumsum(reg_kurtosis_shape) / (1:length(sample_numbers)),
  type = 'l',
  xlab = "Sample Size",
  ylab = "c - Kurtosis",
  col = 'blue',
  lwd = 2
)
plot(
  reg_sample_size,
  cumsum(reg_kurtosis_scale) / (1:length(sample_numbers)),
  type = 'l',
  xlab = "Sample Size",
  ylab = "b - Kurtosis",
  col = 'red',
  lwd = 2
)

par(mfrow = c(2, 2))
plot(
  reg_sample_size,
  mean_shape,
  type = 'l',
  col = 'red',
  lwd = 2,
  xlab = "Sample Size",
  ylab = "c - Shape"
)
plot(
  reg_sample_size,
  mean_scale,
  type = 'l',
  col = 'blue',
  lwd = 2,
  xlab = "Sample Size",
  ylab = "b - Scale"
)
plot(
  reg_sample_size,
  reg_bias_shape,
  type = 'l',
  col = 'red',
  lwd = 2,
  xlab = "Sample Size",
  ylab = "c - Relative Bias"
)
plot(
  reg_sample_size,
  reg_bias_scale,
  type = 'l',
  col = 'blue',
  lwd = 2,
  xlab = "Sample Size",
  ylab = "b - Relative Bias"
)

par(mfrow = c(2, 2))
truehist(
  reg_scale,
  prob = TRUE,
  xlab = "Scale",
  nbins = "fd",
  col = 0
)
curve(
  dnorm(x, mean = mean(reg_scale), sd = sd(reg_scale)),
  from = min(reg_scale) - sd(reg_scale),
  to = max(reg_scale) + sd(reg_scale),
  n = 100,
  add = TRUE,
  type = 'l',
  col = 'blue',
  lwd = 2
)
qqPlot(reg_scale, ylab = "Scale")

truehist(reg_shape,
         prob = T,
         xlab = "Shape",
         col = 0)
curve(
  dnorm(x, mean = mean(reg_shape), sd = sd(reg_shape)),
  from = min(reg_shape) - sd(reg_shape),
  to = max(reg_shape) + sd(reg_shape),
  n = 100,
  add = TRUE,
  type = 'l',
  col = 'blue',
  lwd = 2
)

qqPlot(reg_shape, ylab = "Shape")
