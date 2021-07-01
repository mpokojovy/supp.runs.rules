# (C) Michael Pokojovy and J. Marcus Jobe (2020)

## !!! Set the working directory !!!
setwd("???")
## !!! Set the working directory !!!

set.seed(1)

## IC.ARLs

IC.ARL.WECO   = 91.75
IC.ARL.Nelson = 64.28

CL.WE = qnorm(1 - 0.5/IC.ARL.WECO)
CL.Ne = qnorm(1 - 0.5/IC.ARL.Nelson)

file.name = "data.stream.plots.steady.state=25.pdf"
grDevices::pdf(file.name, width = 10, height = 12)

par(mfcol = c(4, 3))
par(mar = c(4.0, 4.0, 2.0, 0.5))

# Linear trend
slope = .01
IC.time = 1:25
OC.time = 1:64
time = c(IC.time, max(IC.time) + OC.time)

ind = 0L
repeat {
  IC.noise = rnorm(25)
  
  if (max(IC.noise) <= min(CL.WE, CL.Ne)) break
}
OC.noise = rnorm(length(OC.time))
for (sigma in c(1.0, 0.1, 0.0)) {
  ind = ind + 1L
  par(mfg = c(1L, ind))
  
  plot(time, c(IC.noise, slope*OC.time + sigma*OC.noise), 
       xlim = c(min(time), max(time)), ylim = c(min(slope*OC.time) - 3, max(slope*OC.time) + 3),
       xlab = bquote("time period"~t), ylab = bquote(x[t]),
       main = c("Linear drift:", paste("slope = ", slope, ", sigma = ", sigma, sep = "")))
  lines(time, c(IC.noise, slope*OC.time + sigma*OC.noise), lty = 2)
  
  abline(v = max(IC.time), lty = 2)
  
  abline(h =  0, lty = 3)
  abline(h = -CL.WE, lty = 3)
  abline(h =  CL.WE, lty = 3)
  abline(h = -CL.Ne, lty = 3)
  abline(h =  CL.Ne, lty = 3)
  abline(h = -3.0, lty = 3)
  abline(h =  3.0, lty = 3)
}

# Cycling
amplitude = 1.0
period = 16L

IC.time = 1:25
OC.time = 1:64
time = c(IC.time, max(IC.time) + OC.time)

ind = 0L
repeat {
  IC.noise = rnorm(25)
  
  if (max(IC.noise) <= min(CL.WE, CL.Ne)) break
}
OC.noise = rnorm(length(OC.time))
for (sigma in c(1.0, 0.1, 0.0)) {
  ind = ind + 1L
  par(mfg = c(2L, ind))
  
  plot(time, c(IC.noise, amplitude*sin(2*pi*OC.time/period) + sigma*OC.noise), 
       xlim = c(min(time), max(time)), ylim = c(-(amplitude + 3), amplitude + 3),
       xlab = bquote("time period"~t), ylab = bquote(x[t]),
       main = c("Cycling:", paste("amplitude = ", amplitude, ", period = ", period, ", sigma = ", sigma, sep = "")))
  lines(time, c(IC.noise, amplitude*sin(2*pi*OC.time/period) + sigma*OC.noise), lty = 2)
  
  abline(v = max(IC.time), lty = 2)
  
  abline(h =  0, lty = 3)
  abline(h = -CL.WE, lty = 3)
  abline(h =  CL.WE, lty = 3)
  abline(h = -CL.Ne, lty = 3)
  abline(h =  CL.Ne, lty = 3)
  abline(h = -3.0, lty = 3)
  abline(h =  3.0, lty = 3)
}

# Seesaw
amplitude = 0.5
time = 1:64
period = 8L

IC.time = 1:25
OC.time = 1:64
time = c(IC.time, max(IC.time) + OC.time)

ind = 0L
repeat {
  IC.noise = rnorm(25)
  
  if (max(IC.noise) <= min(CL.WE, CL.Ne)) break
}
OC.noise = rnorm(length(OC.time))
for (sigma in c(1.0, 0.1, 0.0)) {
  ind = ind + 1L
  par(mfg = c(3L, ind))
  
  mu = function(t) -amplitude*(-1)^((t - 1) %% period < period/2)
  
  plot(time, c(IC.noise, mu(OC.time) + sigma*OC.noise), 
       xlim = c(min(time), max(time)), ylim = c(-(amplitude + 3), amplitude + 3),
       xlab = bquote("time period"~t), ylab = bquote(x[t]),
       main = c("Seesaw:", paste("amplitude = ", amplitude, ", period = ", period, ", sigma = ", sigma, sep = "")))
  lines(time, c(IC.noise, mu(OC.time) + sigma*OC.noise), lty = 2)
  
  abline(v = max(IC.time), lty = 2)
  
  abline(h =  0, lty = 3)
  abline(h = -CL.WE, lty = 3)
  abline(h =  CL.WE, lty = 3)
  abline(h = -CL.Ne, lty = 3)
  abline(h =  CL.Ne, lty = 3)
  abline(h = -3.0, lty = 3)
  abline(h =  3.0, lty = 3)
}

# Sustained shift
shift = 0.5

IC.time = 1:25
OC.time = 1:64
time = c(IC.time, max(IC.time) + OC.time)

ind = 0L
repeat {
  IC.noise = rnorm(25)
  
  if (max(IC.noise) <= min(CL.WE, CL.Ne)) break
}
OC.noise = rnorm(length(OC.time))
for (sigma in c(1.0, 0.1, 0.0)) {
  ind = ind + 1L
  par(mfg = c(4L, ind))
  
  plot(time, c(IC.noise, shift + sigma*OC.noise), 
       xlim = c(min(time), max(time)), ylim = c(-3, shift + 3),
       xlab = bquote("time period"~t), ylab = bquote(x[t]),
       main = c("Sustained shift:", paste("shift = ", shift, ", sigma = ", sigma, sep = "")))
  lines(time, c(IC.noise, shift + sigma*OC.noise), lty = 2)
  
  abline(v = max(IC.time), lty = 2)
  
  abline(h =  0, lty = 3)
  abline(h = -CL.WE, lty = 3)
  abline(h =  CL.WE, lty = 3)
  abline(h = -CL.Ne, lty = 3)
  abline(h =  CL.Ne, lty = 3)
  abline(h = -3.0, lty = 3)
  abline(h =  3.0, lty = 3)
}

grDevices::dev.off()

try(dev.off(dev.list()["RStudioGD"]), silent = TRUE)
try(dev.off(), silent = TRUE)