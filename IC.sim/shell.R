# (C) Michael Pokojovy and J. Marcus Jobe (2020)

#install.packages("spc")
#install.packages("doParallel")

## !!! Set the working directory !!!
setwd("???")
## !!! Set the working directory !!!

source("../charts.R")
source("simulation.R")

library("spc")
library("doParallel")

seed = 1L
set.seed(seed)

parallel.flag = TRUE

if (parallel.flag) {
  HPC <- makeCluster(detectCores())
  registerDoParallel(HPC)
}

source("run.Nelson.R")
source("run.other.R")
source("plot.results.R")

if (parallel.flag)
    stopCluster(HPC)

save(file = "IC.ARLs.RData", IC.ARL.WECO, IC.ARL.Nelson, IC.stdRL.WECO, IC.stdRL.Nelson)
save(file = "RL.arrays.RData", WECO.RLs, Nelson.RLs, 
     CUSUM.k.0.25.vs.WECO, CUSUM.k.0.50.vs.WECO, CUSUM.k.1.00.vs.WECO,
     CUSUM.k.0.25.vs.Nelson, CUSUM.k.0.50.vs.Nelson, CUSUM.k.1.00.vs.Nelson)