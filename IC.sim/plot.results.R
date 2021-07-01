# (C) Michael Pokojovy and J. Marcus Jobe (2020)

IC.ARL.WECO = 91.75 # According to Champ & Woodall (1987)

for (steady.state.size in c(0L, 25L))
for (pmf.flag in c(FALSE, TRUE)) {
  file.name = if (pmf.flag) 
                paste("RL.pmfs.steady.state.size=", steady.state.size, ".pdf", sep = "") 
              else
                paste("RL.cdfs.steady.state.size=", steady.state.size, ".pdf", sep = "") 
  
  grDevices::pdf(file.name, width = 10, height = 4)

  par(mfcol = c(1, 2))
  par(mar = c(4.0, 4.0, 2.0, 0.5))

  # WECO
  par(mfg = c(1L, 1L))

  p.max = max(1/IC.ARL.WECO, 1/IC.ARL.Nelson)

  RL.grid = seq(from = 1, to = ceiling(5*max(IC.ARL.WECO, IC.ARL.Nelson)), by = 1)

  p = 1/IC.ARL.WECO

  if (pmf.flag) {
    plot(RL.grid, p*(1 - p)^(RL.grid - 1L), xlim = c(1, 5*max(IC.ARL.WECO, IC.ARL.Nelson)), ylim = c(0, p.max*1.2), lty = 2L, lwd = 2, pch = NA, col = 2L,
         xlab = "Run length", ylab = "Probability density function", main = paste("IC ARL = ", zapsmall(IC.ARL.WECO, 4), sep = ""))
    lines(RL.grid, p*(1 - p)^(RL.grid - 1L), lty = 2L, lwd = 2, col = 2L)
    
    RLs = WECO.RLs; RLs = RLs[which(RLs > steady.state.size)] - steady.state.size
    f = density(c(RLs, -(RLs - 1)), from = 1, to = ceiling(5*IC.ARL.WECO))
    lines(f$x, 2.0*f$y, xlab = "Run length", xlim = c(1, 5*max(IC.ARL.WECO, IC.ARL.Nelson)), ylim = c(0, p.max*1.2), lty = 1L, lwd = 2L, col = 1L,
         ylab = "Probability density function", main = paste("IC ARL = ", zapsmall(IC.ARL.WECO, 4), sep = ""))
    
    RLs = CUSUM.k.0.25.vs.WECO; RLs = RLs[which(RLs > steady.state.size)] - steady.state.size
    f = density(c(RLs, -(RLs - 1)), from = 1, to = ceiling(5*IC.ARL.WECO))
    lines(f$x, 2.0*f$y, xlim = c(1, 5*IC.ARL.WECO), lty = 3L, lwd = 2, col = 3L)
    
    RLs = CUSUM.k.0.50.vs.WECO; RLs = RLs[which(RLs > steady.state.size)] - steady.state.size
    f = density(c(RLs, -(RLs - 1)), from = 1, to = ceiling(5*IC.ARL.WECO))
    lines(f$x, 2.0*f$y, xlim = c(1, 5*IC.ARL.WECO), lty = 4L, lwd = 2, col = 4L)
    
    RLs = CUSUM.k.1.00.vs.WECO; RLs = RLs[which(RLs > steady.state.size)] - steady.state.size
    f = density(c(RLs, -(RLs - 1)), from = 1, to = ceiling(5*IC.ARL.WECO))
    lines(f$x, 2.0*f$y, xlim = c(1, 5*IC.ARL.WECO), lty = 5L, lwd = 2, col = 5L)

    legend("topright", legend = c("WE", "Shewhart", "CUSUM(k=0.25)", "CUSUM(k=0.50)", "CUSUM(k=1.00)"),
           lty = 1:5, lwd = 2, col = 1:5)
  } else {
    plot(RL.grid, cumsum(p*(1 - p)^(RL.grid - 1L)), xlab = "Run length", xlim = c(1, 5*max(IC.ARL.WECO, IC.ARL.Nelson)),
         ylim = c(0, 1), lty = 1L, pch = NA, lwd = 2, col = 2L,
         ylab = "Cumulative distribution function", main = paste("IC ARL = ", zapsmall(IC.ARL.WECO, 4), sep = ""))
    lines(RL.grid, cumsum(p*(1 - p)^(RL.grid - 1L)), lty = 1L, lwd = 2, col = 2L)
    
    RLs = WECO.RLs; RLs = RLs[which(RLs > steady.state.size)] - steady.state.size
    lines(RL.grid, ecdf(RLs)(RL.grid), lty = 2L, lwd = 2, col = 1L)
    
    RLs = CUSUM.k.0.25.vs.WECO; RLs = RLs[which(RLs > steady.state.size)] - steady.state.size
    lines(RL.grid, ecdf(RLs)(RL.grid), xlim = c(1, 5*IC.ARL.WECO), lty = 3L, lwd = 2, col = 3L)
    
    RLs = CUSUM.k.0.50.vs.WECO; RLs = RLs[which(RLs > steady.state.size)] - steady.state.size
    lines(RL.grid, ecdf(RLs)(RL.grid), xlim = c(1, 5*IC.ARL.WECO), lty = 4L, lwd = 2, col = 4L)
    
    RLs = CUSUM.k.1.00.vs.WECO; RLs = RLs[which(RLs > steady.state.size)] - steady.state.size
    lines(RL.grid, ecdf(RLs)(RL.grid), xlim = c(1, 5*IC.ARL.WECO), lty = 5L, lwd = 2, col = 5L)

    legend("bottomright", legend = c("WE", "Shewhart", "CUSUM(k=0.25)", "CUSUM(k=0.50)", "CUSUM(k=1.00)"),
           lty = 1:5, lwd = 2, col = 1:5)
  }

  # Nelson
  par(mfg = c(1L, 2L))

  RL.grid = seq(from = 1, to = ceiling(5*max(IC.ARL.WECO, IC.ARL.Nelson)), by = 1)

  p = 1/IC.ARL.Nelson

  if (pmf.flag) {
    plot(RL.grid, p*(1 - p)^(RL.grid - 1L), xlab = "Run length", xlim = c(1, 5*max(IC.ARL.WECO, IC.ARL.Nelson)), ylim = c(0, p.max*1.2), lty = 1L, lwd = 2, pch = NA, col = 1L,
         ylab = "Probability density function", main = paste("IC ARL = ", zapsmall(IC.ARL.Nelson, 4), sep = ""))
    lines(RL.grid, p*(1 - p)^(RL.grid - 1L), lty = 2L, lwd = 2, col = 2L)
    
    RLs = Nelson.RLs; RLs = RLs[which(RLs > steady.state.size)] - steady.state.size
    f = density(c(RLs, -(RLs - 1)), from = 1, to = ceiling(5*IC.ARL.Nelson))
    lines(f$x, 2.0*f$y, xlim = c(1, 5*max(IC.ARL.WECO, IC.ARL.Nelson)), lty = 1L, lwd = 2, col = 1L)
    
    RLs = CUSUM.k.0.25.vs.Nelson; RLs = RLs[which(RLs > steady.state.size)] - steady.state.size
    f = density(c(RLs, -(RLs - 1)), from = 1, to = ceiling(5*IC.ARL.Nelson))
    lines(f$x, 2.0*f$y, xlim = c(1, 5*IC.ARL.Nelson), lty = 3L, lwd = 2, col = 3L)
    
    RLs = CUSUM.k.0.50.vs.Nelson; RLs = RLs[which(RLs > steady.state.size)] - steady.state.size
    f = density(c(RLs, -(RLs - 1)), from = 1, to = ceiling(5*IC.ARL.Nelson))
    lines(f$x, 2.0*f$y, xlim = c(1, 5*IC.ARL.Nelson), lty = 4L, lwd = 2, col = 4L)
    
    RLs = CUSUM.k.1.00.vs.Nelson; RLs = RLs[which(RLs > steady.state.size)] - steady.state.size
    f = density(c(RLs, -(RLs - 1)), from = 1, to = ceiling(5*IC.ARL.Nelson))
    lines(f$x, 2.0*f$y, xlim = c(1, 5*IC.ARL.Nelson), lty = 5L, lwd = 2, col = 5L)

    legend("topright", legend = c("NE", "Shewhart", "CUSUM(k=0.25)", "CUSUM(k=0.50)", "CUSUM(k=1.00)"),
           lty = 1:5, lwd = 2, col = 1:5)
  } else {
    plot(RL.grid, cumsum(p*(1 - p)^(RL.grid - 1L)), xlab = "Run length", xlim = c(1, 5*max(IC.ARL.WECO, IC.ARL.Nelson)),
         ylim = c(0, 1), lty = 1L, pch = NA, lwd = 2, col = 2L,
         ylab = "Cumulative distribution function", main = paste("IC ARL = ", zapsmall(IC.ARL.Nelson, 4), sep = ""))
    lines(RL.grid, cumsum(p*(1 - p)^(RL.grid - 1L)), lty = 1L, lwd = 2, col = 2L)
    
    RLs = Nelson.RLs; RLs = RLs[which(RLs > steady.state.size)] - steady.state.size
    lines(RL.grid, ecdf(RLs)(RL.grid), lty = 2L, lwd = 2, col = 1L)
    
    RLs = CUSUM.k.0.25.vs.Nelson; RLs = RLs[which(RLs > steady.state.size)] - steady.state.size
    lines(RL.grid, ecdf(RLs)(RL.grid), xlim = c(1, 5*IC.ARL.WECO), lty = 3L, lwd = 2, col = 3L)
    
    RLs = CUSUM.k.0.50.vs.Nelson; RLs = RLs[which(RLs > steady.state.size)] - steady.state.size
    lines(RL.grid, ecdf(RLs)(RL.grid), xlim = c(1, 5*IC.ARL.WECO), lty = 4L, lwd = 2, col = 4L)
    
    RLs = CUSUM.k.1.00.vs.Nelson; RLs = RLs[which(RLs > steady.state.size)] - steady.state.size
    lines(RL.grid, ecdf(RLs)(RL.grid), xlim = c(1, 5*IC.ARL.WECO), lty = 5L, lwd = 2, col = 5L)

    legend("bottomright", legend = c("NE", "Shewhart", "CUSUM(k=0.25)", "CUSUM(k=0.50)", "CUSUM(k=1.00)"),
           lty = 1:5, lwd = 2, col = 1:5)
  }

  grDevices::dev.off()
}

try(dev.off(dev.list()["RStudioGD"]), silent = TRUE)
try(dev.off(), silent = TRUE)