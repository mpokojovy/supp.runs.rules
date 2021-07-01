# (C) Michael Pokojovy and J. Marcus Jobe (2020)

#install.packages("spc")
#install.packages("doParallel")

## !!! Set the working directory !!!
setwd("???")
## !!! Set the working directory !!!

source("../charts.R")
source("simulation.R")
source("auxil.R")

library("spc")
library("doParallel")

##

seed = 1L
set.seed(seed)

parallel.flag = TRUE

if (parallel.flag) {
  HPC <- makeCluster(detectCores())
  registerDoParallel(HPC)
}

## simulations settings
nrep = 1000L # small-scale simulation
nmax = 5000L

## settings: sigma_IC = 1 by default
sigma.new = 1.0 # sigma_OC

steady.state.sizes = c(0L, 25L)

ks = c(0.25, 0.50, 1.00)

slopes      = c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0)
amplitudes  = c(0.25, 0.50, 0.75, 1.00, 1.50, 2.0, 2.75, 3.0, 3.5)
periods     = c(4L, 8L, 16L, 32L, 64L)
shifts      = c(0.25, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00, 3.5)

n.scenarios = (2L + 2L + 2*length(ks))*length(steady.state.sizes)*(length(slopes) + 2L*length(amplitudes)*length(periods) + length(shifts))

results.df = data.frame(matrix(0.0, ncol = 16L, nrow = n.scenarios))
colnames(results.df) = c("chart.index", "k", "IC.ARL", "steady.state.size", "ARL", "VRL", "std.RL", "ARL0", "VRL0", "std.RL0",
                         "surv.rate", "overrun.rate", "sigma.new", "shift.type", "par.1", "par.2")

## plotting settings
plot.capped = TRUE

## IC.ARLs

IC.ARL.WECO   = 91.75
IC.ARL.Nelson = 64.28

cat("IC.ARL.WECO = ", IC.ARL.WECO, "\n", sep = "")
cat("IC.ARL.Nelson = ", IC.ARL.Nelson, "\n", sep = "")

## Run

ind.scenario = 0L

ptm = proc.time()

for (steady.state.size in steady.state.sizes)
for (chart.index in 0:3) {
  .ks      = if (chart.index == 3) ks else NA
  .IC.ARLs = if (chart.index == 0) IC.ARL.WECO else {if (chart.index == 1) IC.ARL.Nelson else c(IC.ARL.WECO, IC.ARL.Nelson)}
  
  for (IC.ARL in .IC.ARLs)
  for (k in .ks) {
    # linear shift
    shift.type = 1L
    
    for (slope in slopes) {
      ind.scenario = ind.scenario + 1L
      
      results.df[ind.scenario, ] = do.simulation(chart.index = chart.index, k = k, IC.ARL = IC.ARL, steady.state.size = steady.state.size, 
                                                 nrep = nrep, nmax = nmax, sigma.new = sigma.new, shift.type = shift.type, par = c(slope))
    }
    
    # cycling
    shift.type = 2L
    
    for (amplitude in amplitudes)
    for (period in periods) {
      ind.scenario = ind.scenario + 1L
      
      
      results.df[ind.scenario, ] = do.simulation(chart.index = chart.index, k = k, IC.ARL = IC.ARL, steady.state.size = steady.state.size, 
                                                 nrep = nrep, nmax = nmax, sigma.new = sigma.new, shift.type = shift.type, par = c(amplitude, period))
    }
    
    # seesaw
    shift.type = 3L
    
    for (amplitude in amplitudes)
    for (period in periods) {
      ind.scenario = ind.scenario + 1L
      
      results.df[ind.scenario, ] = do.simulation(chart.index = chart.index, k = k, IC.ARL = IC.ARL, steady.state.size = steady.state.size, 
                                                 nrep = nrep, nmax = nmax, sigma.new = sigma.new, shift.type = shift.type, par = c(amplitude, period))
    }
    
    # shift
    shift.type = 4L
    
    for (shift in shifts) {
      ind.scenario = ind.scenario + 1L
      
      results.df[ind.scenario, ] = do.simulation(chart.index = chart.index, k = k, IC.ARL = IC.ARL, steady.state.size = steady.state.size, 
                                                   nrep = nrep, nmax = nmax, sigma.new = sigma.new, shift.type = shift.type, par = c(shift))
    }
  }
}

print(proc.time() - ptm)

if (parallel.flag)
  stopCluster(HPC)

save(file = paste("results.sigma.new=", sigma.new, ".seed=", seed, ".RData", sep = ""), 
     results.df, nrep, nmax, sigma.new, steady.state.sizes, ks, slopes, amplitudes, periods, shifts)

table.index = 0L # table index offset

results.as.text()

if (!dir.exists("fig"))
  dir.create("fig/")  

plot.results()