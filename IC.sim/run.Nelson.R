# (C) Michael Pokojovy and J. Marcus Jobe (2020)

## simulations settings
nrep    = 100E6L
nchunks = 10000L
nmax    = 1E6L

Nelson.RLs = rep(0L, nrep)

ptm = proc.time()

for (chunk in 1:nchunks) {
  cat("Nelson: chunk = ", chunk, "\n", sep = "")
  I = ((chunk - 1)*(nrep/nchunks) + 1):(chunk*(nrep/nchunks))
  
  ## Nelson
  Nelson.RLs[I] = do.simulation(chart.index = 1L, k = NA, IC.ARL = NA, nrep = nrep/nchunks, nmax = nmax, shift.type = 0L, par = NULL)
}

IC.ARL.Nelson   = mean(Nelson.RLs)
IC.stdRL.Nelson = sd(Nelson.RLs)

print(proc.time() - ptm)

cat("IC ARL Nelson = ", IC.ARL.Nelson, "\n", sep = "")
cat("IC std RL Nelson = ", IC.stdRL.Nelson, "\n", sep = "")
cat("IC std error ARL Nelson = ", IC.stdRL.Nelson/sqrt(nrep), "\n", sep = "")