# (C) Michael Pokojovy and J. Marcus Jobe (2020)

## simulations settings
nrep    = 10E6L
nchunks = 1000L
nmax    = 1E6L

WECO.RLs = rep(0L, nrep)

ptm = proc.time()

for (chunk in 1:nchunks) {
  cat("WECO: chunk = ", chunk, "\n", sep = "")
  I = ((chunk - 1)*(nrep/nchunks) + 1):(chunk*(nrep/nchunks))
  
  ## Western Electric
  WECO.RLs[I]   = do.simulation(chart.index = 0L, k = NA, IC.ARL = NA, nrep = nrep/nchunks, nmax = nmax, shift.type = 0L, par = NULL)
}

print(proc.time() - ptm)

IC.ARL.WECO     = mean(WECO.RLs)
IC.stdRL.WECO   = sd(WECO.RLs)

cat("IC ARL WECO = ", IC.ARL.WECO, "\n", sep = "")
cat("IC std RL WECO = ", IC.stdRL.WECO, "\n", sep = "")
cat("IC std error ARL WECO = ", IC.stdRL.WECO/sqrt(nrep), "\n", sep = "")

CUSUM.k.0.25.vs.WECO   = rep(0L, nrep)
CUSUM.k.0.50.vs.WECO   = rep(0L, nrep)
CUSUM.k.1.00.vs.WECO   = rep(0L, nrep)
CUSUM.k.0.25.vs.Nelson = rep(0L, nrep)
CUSUM.k.0.50.vs.Nelson = rep(0L, nrep)
CUSUM.k.1.00.vs.Nelson = rep(0L, nrep)

for (chunk in 1:nchunks) {
  cat("CUSUM: chunk = ", chunk, "\n", sep = "")
  I = ((chunk - 1)*(nrep/nchunks) + 1):(chunk*(nrep/nchunks))
  
  ## CUSUM
  CUSUM.k.0.25.vs.WECO[I] = do.simulation(chart.index = 3L, k = 0.25, IC.ARL = IC.ARL.WECO, nrep = nrep/nchunks, nmax = nmax, shift.type = 0L, par = NULL)
  CUSUM.k.0.50.vs.WECO[I] = do.simulation(chart.index = 3L, k = 0.50, IC.ARL = IC.ARL.WECO, nrep = nrep/nchunks, nmax = nmax, shift.type = 0L, par = NULL)
  CUSUM.k.1.00.vs.WECO[I] = do.simulation(chart.index = 3L, k = 1.00, IC.ARL = IC.ARL.WECO, nrep = nrep/nchunks, nmax = nmax, shift.type = 0L, par = NULL)
  
  CUSUM.k.0.25.vs.Nelson[I] = do.simulation(chart.index = 3L, k = 0.25, IC.ARL = IC.ARL.Nelson, nrep = nrep/nchunks, nmax = nmax, shift.type = 0L, par = NULL)
  CUSUM.k.0.50.vs.Nelson[I] = do.simulation(chart.index = 3L, k = 0.50, IC.ARL = IC.ARL.Nelson, nrep = nrep/nchunks, nmax = nmax, shift.type = 0L, par = NULL)
  CUSUM.k.1.00.vs.Nelson[I] = do.simulation(chart.index = 3L, k = 1.00, IC.ARL = IC.ARL.Nelson, nrep = nrep/nchunks, nmax = nmax, shift.type = 0L, par = NULL)
}

print(proc.time() - ptm)