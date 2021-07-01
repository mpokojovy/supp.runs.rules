# (C) Michael Pokojovy and J. Marcus Jobe (2020)

zapsmall3 <- function(x) {
  if (length(x) == 0) {
    return("!!!")
  }

  if (x == 0)
    return("0")

  if ((abs(x) > 1E-3) && (abs(x) < 1E3)) {
    return(sprintf("%4.3f", x))
  } else {
    n = log10(abs(x))
    if (n >= 0) {
      n = ceiling(n)
      c = round(x*10^(-n), digits = 0)
      if (c == 0) n = n - 1
    } else {
      n = floor(n)
      c = round(x*10^(-n), digits = 0)
      if (c == 0) n = n + 1
    }

    c = round(x*10^(-n), digits = 0)

    if (abs(c) < 10) {
      return(paste(c, "\\text{E", n, "}", sep = ""))
    } else if (n > 0) {
      c = c/10
      n = n + 1
      return(paste(c, "\\text{E", n, "}", sep = ""))
    } else {
      c = c/10
      n = n - 1
      return(paste(c, "\\text{E", n, "}", sep = ""))
    }
  }
}

get.ARL <- function(chart.index, k, IC.ARL, steady.state.size, sigma.new, shift.type, par, stat = 0) {
  if (shift.type == 0) {
    if (chart.index == 0L) {
      return(IC.ARL.WECO)
    } else if (chart.index == 1L) {
      return(IC.ARL.Nelson)
    } else {
      return(IC.ARL)
    }
  }
  
  I = which((results.df$chart.index == chart.index) & (results.df$steady.state.size == steady.state.size) & 
              (results.df$sigma.new == sigma.new) & (results.df$shift.type == shift.type))
  
  if (chart.index > 1L) {
    I = I[which(results.df[I, ]$IC.ARL == IC.ARL)]
    
    if (chart.index > 2L) {
      I = I[which(results.df[I, ]$k == k)]
    }
  }
  
  if (shift.type > 0L) {
    I = I[which(results.df[I, ]$par.1 == par[1])]
    
    if ((shift.type > 1L) && (shift.type < 4L)) {
      I = I[which(results.df[I, ]$par.2 == par[2])]
    }
  }
  
  if (length(I) == 0) {
    return(NA)
  } else {
    if (stat == 0)
      return(results.df$ARL[I[1]])
    else
      return(results.df$std.RL[I[1]]/sqrt(nrep*results.df$surv.rate[I[1]]))
  }
}

get.ARL.interval <- function(chart.index, k, IC.ARL, steady.state.size, sigma.new, shift.type, par, stat = 0) {
  if (shift.type == 0) {
    if (chart.index == 0L) {
      return(IC.ARL.WECO)
    } else if (chart.index == 1L) {
      return(IC.ARL.Nelson)
    } else {
      return(IC.ARL)
    }
  }
  
  I = which((results.df$chart.index == chart.index) & (results.df$steady.state.size == steady.state.size) & 
              (results.df$sigma.new == sigma.new) & (results.df$shift.type == shift.type))
  
  if (chart.index > 1L) {
    I = I[which(results.df[I, ]$IC.ARL == IC.ARL)]
    
    if (chart.index > 2L) {
      I = I[which(results.df[I, ]$k == k)]
    }
  }
  
  if (shift.type > 0L) {
    I = I[which(results.df[I, ]$par.1 == par[1])]
    
    if ((shift.type > 1L) && (shift.type < 4L)) {
      I = I[which(results.df[I, ]$par.2 == par[2])]
    }
  }
  
  if (length(I) == 0) {
    return("")
  } else {
    n = nrep*results.df$surv.rate[I[1]]
    
    if ((results.df$steady.state.size[I[1]] == 0) && (results.df$sigma.new[I[1]] == 0)) {
      ARL = results.df$ARL[I[1]]
      ARL = max(ARL, 1.0)
      
      if (results.df$overrun.rate[I[1]] == 0) {
        return(paste(zapsmall3(ARL), sep = ""))
      } else {
        return(paste("\\geq", zapsmall3(ARL), sep = ""))
      }
    } else {
      if (results.df$overrun.rate[I[1]] == 0) {
        ARL     = results.df$ARL[I[1]]
        ARL     = max(ARL, 1.0)
        std.ARL = results.df$std.RL[I[1]]/sqrt(n)
        
        return(paste(zapsmall3(ARL), "\\pm 3(", zapsmall3(std.ARL), ")", sep = ""))
      } else if (results.df$ARL[I[1]] >= 1E3 - 1E-6) {
        return("\\geq 1000")
      } else {
        # p = results.df$surv.rate[I[1]]*results.df$overrun.rate[I[1]]
        # 
        # ARL0     = results.df$ARL0[I[1]]
        # std.ARL0 = results.df$std.RL0[I[1]]/sqrt(n)
        
        ARL     = results.df$ARL[I[1]]
        ARL     = max(ARL, 1.0)
        std.ARL = results.df$std.RL[I[1]]/sqrt(n)
        
        cat("ARL =", ARL, "std.ARL =", std.ARL, "\n")
        
        return(paste("\\geq", zapsmall3(ARL), "- 3(", zapsmall3(std.ARL), ")", sep = ""))
      }
    }
  }
}

get.ARL.comp <- function(chart.index1, chart.index2, k, IC.ARL, steady.state.size, sigma.new, shift.type, par) {
  return(if (is.na(chart.index2)) get.ARL(chart.index1, k, IC.ARL, steady.state.size, sigma.new, shift.type, par)
         else get.ARL(chart.index1, k, IC.ARL, steady.state.size, sigma.new, shift.type, par)/
           get.ARL(chart.index2, k, IC.ARL, steady.state.size, sigma.new, shift.type, par))  
}

results.as.text <- function() {
  file = paste("comparisons.sigma.new=", sigma.new, ".txt", sep = "")
  cat("", file = file, append = FALSE)
  
  #cat("\\section{ARL Confidence Interval Tables for $\\sigma_{\\mathrm{OC}} = ", sigma.new, "$}\n", file = file, append = TRUE)
  
  for (steady.state.size in steady.state.sizes) {
    cat("\\subsection{Steady-State Size ", steady.state.size, "}\n", file = file, append = TRUE)
    
    ## Linear
    cat("\\begin{longtable}{c|", paste(rep("c", (2 + length(ks))), collapse = ""), "}\n", file = file, append = TRUE)
    
    for (IC.ARL in c(IC.ARL.WECO, IC.ARL.Nelson)) {
      cat("\t\t\\hline\n", file = file, append = TRUE)
      
      cat("\t\t & \\multicolumn{", 2 + length(ks), "}{c}{$\\mathrm{IC\\,ARL} = ", IC.ARL, "$} \\\\ \n", file = file, append = TRUE)
      
      row = paste("\t\tSlope $\\Delta$ &", if (IC.ARL == IC.ARL.WECO) "WE" else "NE", " & Shewhart $X$ & ",
                  paste(sapply(ks, function(k) paste(" CUSUM($k=", k, "$) ", sep = "")), collapse = "&"))
      
      cat(row, "\\\\ \n", file = file, append = TRUE)
      cat("\t\t\\hline \n", file = file, append = TRUE)
      
      for (slope in slopes) {
        row = ""
        
        interval = get.ARL.interval(if (IC.ARL == IC.ARL.WECO) 0L else 1L, k = NA, IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 1, c(slope, NA), stat = 0)
        
        row = paste(row, "$", interval, "$ & ", sep = "")
        
        interval = get.ARL.interval(2L, k = NA, IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 1, c(slope, NA), stat = 0)
        
        row = paste(row, "$", interval, "$ ", sep = "")
        
        for (k in ks) {
          interval = get.ARL.interval(3L, k, IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 1, c(slope, NA), stat = 0)
          
          row = paste(row, "& $", interval, "$ ", sep = "")
        }
        
        cat("\t\t", slope, " & ", row, "\\\\ \n", file = file, append = TRUE)
      }
      
      cat("\t\t\\hline \n", file = file, append = TRUE)
    }
    
    table.index = table.index + 1L
    caption = paste("ARL confidence intervals, linear drift, $\\sigma_{\\mathrm{OC}} = ", sigma.new, "$, ", 
                    "steady state size ", steady.state.size, ", cf.~Supplement Figure \\ref{FIGURE:SUPPLEMENT_", sprintf("%02d", table.index + 2L), "}. ", 
                    "\\label{TABLE:SUPPLEMENT_", table.index, "}", sep = "")
    
    cat("\\caption{", caption, "}", file = file, append = TRUE)
    cat("\\end{longtable}\n", file = file, append = TRUE)
    cat("\n\n\\newpage\n\n", file = file, append = TRUE)
    
    ## Sustained
    cat("\\begin{longtable}{c|", paste(rep("c", (2 + length(ks))), collapse = ""), "}\n", file = file, append = TRUE)
    
    for (IC.ARL in c(IC.ARL.WECO, IC.ARL.Nelson)) {
      cat("\t\t\\hline\n", file = file, append = TRUE)
      
      cat("\t\t & \\multicolumn{", 2 + length(ks), "}{c}{$\\mathrm{IC\\,ARL} = ", IC.ARL, "$} \\\\ \n", file = file, append = TRUE)
      
      row = paste("\t\tShift $\\mu_{1}$ &", if (IC.ARL == IC.ARL.WECO) "WE" else "NE", " & Shewhart $X$ & ",
                  paste(sapply(ks, function(k) paste(" CUSUM($k=", k, "$) ", sep = "")), collapse = "&"))
      
      cat(row, "\\\\ \n", file = file, append = TRUE)
      cat("\t\t\\hline \n", file = file, append = TRUE)
      
      for (shift in shifts) {
        row = ""
        
        interval = get.ARL.interval(if (IC.ARL == IC.ARL.WECO) 0L else 1L, k = NA, IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 4, c(shift, NA), stat = 0)

        row = paste(row, "$", interval, "$ & ", sep = "")
        
        interval = get.ARL.interval(2L, k = NA, IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 4, c(shift, NA), stat = 0)
        
        row = paste(row, "$", interval, "$ ", sep = "")
        
        for (k in ks) {
          interval = get.ARL.interval(3L, k, IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 4, c(shift, NA), stat = 0)
          
          row = paste(row, "& $", interval, "$ ", sep = "")
        }
        
        cat("\t\t", shift, " & ", row, "\\\\ \n", file = file, append = TRUE)
      }
      
      cat("\t\t\\hline \n", file = file, append = TRUE)
    }
    
    table.index = table.index + 1L
    caption = paste("ARL confidence intervals, sustained shift, $\\sigma_{\\mathrm{OC}} = ", sigma.new, "$, ", 
                    "steady state size ", steady.state.size, ", cf.~Supplement Figure \\ref{FIGURE:SUPPLEMENT_", sprintf("%02d", table.index + 2L), "}. ", 
                    "\\label{TABLE:SUPPLEMENT_", table.index, "}", sep = "")
    
    cat("\\caption{", caption, "}", file = file, append = TRUE)
    cat("\\end{longtable}\n", file = file, append = TRUE)
    cat("\n\n\\newpage\n\n", file = file, append = TRUE)
    
    ## Cycling
    cat("\\begin{longtable}{cc|", paste(rep("c", (2 + length(ks))), collapse = ""), "}\n", file = file, append = TRUE)
    for (IC.ARL in c(IC.ARL.WECO, IC.ARL.Nelson)) {
      cat("\t\t\\hline\n", file = file, append = TRUE)
      cat("\t\t & & \\multicolumn{", 2 + length(ks), "}{c}{$\\mathrm{IC\\,ARL} = ", IC.ARL, "$} \\\\ \n", file = file, append = TRUE)

      row = paste("\t\tperiod $\\omega$ & ampl. $A$ & ", if (IC.ARL == IC.ARL.WECO) "WE" else "NE", " & Shewhart $X$ & ",
                  paste(sapply(ks, function(k) paste(" CUSUM($k=", k, "$) ", sep = "")), collapse = "&"))
      
      cat(row, "\\\\ \n", file = file, append = TRUE)
      cat("\t\t\\hline \n", file = file, append = TRUE)
      
      for (period in periods)
        for (amplitude in amplitudes) {
          row = ""
          
          interval = get.ARL.interval(if (IC.ARL == IC.ARL.WECO) 0L else 1L, k = NA, IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 2, c(amplitude, period), stat = 0)

          row = paste(row, "$", interval, "$ & ", sep = "")
          
          interval = get.ARL.interval(2L, k = NA, IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 2, c(amplitude, period), stat = 0)
          
          row = paste(row, "$", interval, "$ ", sep = "")
          
          for (k in ks) {
            interval = get.ARL.interval(3L, k, IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 2, c(amplitude, period), stat = 0)
            
            row = paste(row, "& $", interval, "$ ", sep = "")
          }
          
          cat("\t\t", period, " & ", amplitude, " & ", row, "\\\\ \n", file = file, append = TRUE)
        }
      
      cat("\t\t\\hline \n", 
          if (IC.ARL == IC.ARL.WECO) 
            "\\\\ \\caption*{(Continued on the next page)} \\\\ \\pagebreak" else "", file = file, append = TRUE)
    }
    
    table.index = table.index + 1L
    caption = paste("ARL confidence intervals, cycling, $\\sigma_{\\mathrm{OC}} = ", sigma.new, "$, ", 
                    "steady state size ", steady.state.size, ", cf.~Supplement Figure \\ref{FIGURE:SUPPLEMENT_", sprintf("%02d", table.index + 2L), "}. ", 
                    "\\label{TABLE:SUPPLEMENT_", table.index, "}", sep = "")
    
    cat("\\caption{", caption, "}", file = file, append = TRUE)
    cat("\\end{longtable}\n", file = file, append = TRUE)
    cat("\n\n\\newpage\n\n", file = file, append = TRUE)
    
    ## Seesaw
    cat("\\begin{longtable}{cc|", paste(rep("c", (2 + length(ks))), collapse = ""), "}\n", file = file, append = TRUE)
    for (IC.ARL in c(IC.ARL.WECO, IC.ARL.Nelson)) {
      cat("\t\t\\hline\n", file = file, append = TRUE)
      cat("\t\t & & \\multicolumn{", 2 + length(ks), "}{c}{$\\mathrm{IC\\,ARL} = ", IC.ARL, "$} \\\\ \n", file = file, append = TRUE)
      
      row = paste("\t\tperiod $\\omega$ & ampl. $A$ & ", if (IC.ARL == IC.ARL.WECO) "WE" else "NE", " & Shewhart $X$ & ",
                  paste(sapply(ks, function(k) paste(" CUSUM($k=", k, "$) ", sep = "")), collapse = "&"))
      
      cat(row, "\\\\ \n", file = file, append = TRUE)
      cat("\t\t\\hline \n", file = file, append = TRUE)
      
      for (period in periods)
        for (amplitude in amplitudes) {
          row = ""
          
          interval = get.ARL.interval(if (IC.ARL == IC.ARL.WECO) 0L else 1L, k = NA, IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 3, c(amplitude, period), stat = 0)

          row = paste(row, "$", interval, "$ & ", sep = "")
          
          interval = get.ARL.interval(2L, k = NA, IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 3, c(amplitude, period), stat = 0)

          row = paste(row, "$", interval, "$ ", sep = "")
          
          for (k in ks) {
            interval = get.ARL.interval(3L, k, IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 3, c(amplitude, period), stat = 0)
            
            row = paste(row, "& $", interval, "$ ", sep = "")
          }
          
          cat("\t\t", period, " & ", amplitude, " & ", row, "\\\\ \n", file = file, append = TRUE)
        }
      
      cat("\t\t\\hline \n", 
          if (IC.ARL == IC.ARL.WECO) 
            "\\\\ \\caption*{(Continued on the next page)} \\\\ \\pagebreak" else "", file = file, append = TRUE)
    }
    
    table.index = table.index + 1L
    caption = paste("ARL confidence intervals, seesaw, $\\sigma_{\\mathrm{OC}} = ", sigma.new, "$, ", 
                    "steady state size ", steady.state.size, ", cf.~Supplement Figure \\ref{FIGURE:SUPPLEMENT_", sprintf("%02d", table.index + 2L), "}. ", 
                    "\\label{TABLE:SUPPLEMENT_", table.index, "}", sep = "")
    
    cat("\\caption{", caption, "}", file = file, append = TRUE)
    cat("\\end{longtable}\n", file = file, append = TRUE)
    cat("\n\n\\newpage\n\n", file = file, append = TRUE)
  }
}

plot.results <- function() {
  for (steady.state.size in steady.state.sizes) {
    ## Linear
    file.name = paste("fig/linear.ss.size=", steady.state.size, ".sigma.new=", sigma.new, ".pdf", sep = "")
    grDevices::pdf(file.name, width = 10, height = 4)
    
    par(mfcol = c(1, 2))
    par(mar = c(4.0, 4.0, 2.0, 0.5))
    
    for (IC.ARL in c(IC.ARL.WECO, IC.ARL.Nelson)) {
      runs.rules.curve = rep(0.0, length(slopes))
      shewhart.curve = rep(0.0, length(slopes))
      cusum.curves = matrix(0.0, nrow = length(ks), ncol = length(slopes))
      
      par(mfg = c(1L, if (IC.ARL == IC.ARL.WECO) 1L else 2L))
      
      for (j in 1:length(slopes)) {
        runs.rules.curve[j] = get.ARL(if (IC.ARL == IC.ARL.WECO) 0L else 1L, k = NA, IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 1L, c(slopes[j], NA))
        shewhart.curve[j]   = get.ARL(2L, k = NA, IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 1L, c(slopes[j], NA))
        
        for (i in 1:length(ks)) {
          cusum.curves[i, j] = get.ARL(3L, k = ks[i], IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 1L, c(slopes[j], NA))
        }
      }
      
      xmin = min(log(slopes))
      xmax = max(log(slopes))
      ymin = 0.0
      ymax = if (plot.capped) IC.ARL*1.2 else max(max(runs.rules.curve, shewhart.curve, cusum.curves)*1.5, IC.ARL*2.0)
      
      plot(log(slopes), runs.rules.curve, xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = bquote(log("slope")), ylab = "ARL",
           xaxt = "n",
           type = "o", lty = 1, pch = 1, col = "black", 
           main = c("Linear trending",
                    (if (IC.ARL == IC.ARL.WECO) "WE vs. Shewhart vs. CUSUM" else "NE vs. Shewhart vs. CUSUM")))
      axis(side = 1, at = log(slopes), labels = sapply(slopes, function(s) paste("log(", s, ")", sep = "")))
      
      lines(log(slopes), shewhart.curve, type = "o", lty = 1, pch = 0, col = "red")
      sapply(1:length(ks), function(i) lines(log(slopes), as.vector(cusum.curves[i, ]), type = "o", lty = i, pch = 2, col = "blue"))
      abline(h = IC.ARL, lty = 2)
      
      legend("topright", legend = c("IC ARL", if (IC.ARL == IC.ARL.WECO) "WE" else "NE", "Shewhart",
                                    sapply(1:length(ks), function(i) paste("CUSUM(k=", ks[i], ")", sep = ""))),
             pch = c(NA, 1, 0, rep(2, length(ks))),
             lty = c(2, 1, 1, 1:length(ks)),
             col = c("black", "black", "red", rep("blue", length(ks))), bg = "white", pt.cex = 1, cex = 0.7)
    }
    
    grDevices::dev.off()
    
    # Cycling
    file.name = paste("fig/cycling.ss.size=", steady.state.size, ".sigma.new=", sigma.new, ".pdf", sep = "")
    grDevices::pdf(file.name, width = 8, height = 14*length(periods)/5)
    
    par(mfcol = c(length(periods), 2))
    par(mar = c(4.0, 4.0, 2.0, 0.5))
    
    for (IC.ARL in c(IC.ARL.WECO, IC.ARL.Nelson))
      for (i.period in 1:length(periods)) {
        par(mfg = c(i.period, if (IC.ARL == IC.ARL.WECO) 1L else 2L))
        
        period = periods[i.period]
        
        runs.rules.curve = rep(0.0, length(amplitudes))
        shewhart.curve = rep(0.0, length(amplitudes))
        cusum.curves = matrix(0.0, nrow = length(ks), ncol = length(amplitudes))
        
        for (j in 1:length(amplitudes)) {
          runs.rules.curve[j] = get.ARL(if (IC.ARL == IC.ARL.WECO) 0L else 1L, k = NA, IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 2L, c(amplitudes[j], period))
          shewhart.curve[j]   = get.ARL(2L, k = NA, IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 2L, c(amplitudes[j], period))
          
          for (i in 1:length(ks)) {
            cusum.curves[i, j] = get.ARL(3L, k = ks[i], IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 2L, c(amplitudes[j], period))
          }
        }
        
        xmin = min(amplitudes)
        xmax = max(amplitudes)
        ymin = 0.0
        ymax = if (plot.capped) IC.ARL*1.2 else max(max(runs.rules.curve, shewhart.curve, cusum.curves)*1.5, IC.ARL*2.0)
        
        plot(amplitudes, runs.rules.curve, xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = "amplitude", ylab = "ARL",
             xaxt = "n",
             type = "o", lty = 1, pch = 1, col = "black",
             main = c(paste("Cycling: period =", period),
                      (if (IC.ARL == IC.ARL.WECO) "WE vs. Shewhart vs. CUSUM" else "NE vs. Shewhart vs. CUSUM")))
        axis(side = 1, at = amplitudes, labels = sapply(amplitudes, function(s) paste(zapsmall(s, 2), sep = "")))
        
        lines(amplitudes, shewhart.curve, type = "o", lty = 1, pch = 0, col = "red")
        sapply(1:length(ks), function(i) lines(amplitudes, as.vector(cusum.curves[i, ]), type = "o", lty = i, pch = 2, col = "blue"))
        abline(h = IC.ARL, lty = 2)
        
        legend("topright", legend = c("IC ARL", if (IC.ARL == IC.ARL.WECO) "WE" else "NE", "Shewhart",
                                      sapply(1:length(ks), function(i) paste("CUSUM(k=", ks[i], ")", sep = ""))),
               pch = c(NA, 1, 0, rep(2, length(ks))),
               lty = c(2, 1, 1, 1:length(ks)),
               col = c("black", "black", "red", rep("blue", length(ks))), bg = "white", pt.cex = 1, cex = 0.7)
      }
    
    grDevices::dev.off()
    
    # Seesaw
    file.name = paste("fig/seesaw.ss.size=", steady.state.size, ".sigma.new=", sigma.new, ".pdf", sep = "")
    grDevices::pdf(file.name, width = 8, height = 14*length(periods)/5)
    
    par(mfcol = c(length(periods), 2))
    par(mar = c(4.0, 4.0, 2.0, 0.5))
    
    for (IC.ARL in c(IC.ARL.WECO, IC.ARL.Nelson))
      for (i.period in 1:length(periods)) {
        par(mfg = c(i.period, if (IC.ARL == IC.ARL.WECO) 1L else 2L))
        
        period = periods[i.period]
        
        runs.rules.curve = rep(0.0, length(amplitudes))
        shewhart.curve = rep(0.0, length(amplitudes))
        cusum.curves = matrix(0.0, nrow = length(ks), ncol = length(amplitudes))
        
        for (j in 1:length(amplitudes)) {
          runs.rules.curve[j] = get.ARL(if (IC.ARL == IC.ARL.WECO) 0L else 1L, k = NA, IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 3L, c(amplitudes[j], period))
          shewhart.curve[j] = get.ARL(2L, k = NA, IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 3L, c(amplitudes[j], period))
          
          for (i in 1:length(ks)) {
            cusum.curves[i, j] = get.ARL(3L, k = ks[i], IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 3L, c(amplitudes[j], period))
          }
        }
        
        xmin = min(amplitudes)
        xmax = max(amplitudes)
        ymin = 0.0
        ymax = if (plot.capped) IC.ARL*1.2 else max(max(runs.rules.curve, shewhart.curve, cusum.curves)*1.5, IC.ARL*2.0)
        
        plot(amplitudes, runs.rules.curve, xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = "amplitude", ylab = "ARL",
             xaxt = "n",
             type = "o", lty = 1, pch = 1, col = "black",
             main = c(paste("Seesaw upset: period =", period),
                      (if (IC.ARL == IC.ARL.WECO) "WE vs. Shewhart vs. CUSUM" else "NE vs. Shewhart vs. CUSUM")))
        axis(side = 1, at = amplitudes, labels = sapply(amplitudes, function(s) paste(s, sep = "")))
        
        lines(amplitudes, shewhart.curve, type = "o", lty = 1, pch = 0, col = "red")
        sapply(1:length(ks), function(i) lines(amplitudes, as.vector(cusum.curves[i, ]), type = "o", lty = i, pch = 2, col = "blue"))
        abline(h = IC.ARL, lty = 2)
        
        legend("topright", legend = c("IC ARL", if (IC.ARL == IC.ARL.WECO) "WE" else "NE", "Shewhart",
                                      sapply(1:length(ks), function(i) paste("CUSUM(k=", ks[i], ")", sep = ""))),
               pch = c(NA, 1, 0, rep(2, length(ks))),
               lty = c(2, 1, 1, 1:length(ks)),
               col = c("black", "black", "red", rep("blue", length(ks))), bg = "white", pt.cex = 1, cex = 0.7)
      }
    
    grDevices::dev.off()
    
    ## Sustained
    file.name = paste("fig/sustained.ss.size=", steady.state.size, ".sigma.new=", sigma.new, ".pdf", sep = "")
    grDevices::pdf(file.name, width = 10, height = 4)
    
    par(mfcol = c(1, 2))
    par(mar = c(4.0, 4.0, 2.0, 0.5))
    
    for (IC.ARL in c(IC.ARL.WECO, IC.ARL.Nelson)) {
      runs.rules.curve = rep(0.0, length(shifts))
      shewhart.curve = rep(0.0, length(shifts))
      cusum.curves = matrix(0.0, nrow = length(ks), ncol = length(shifts))
      
      par(mfg = c(1L, if (IC.ARL == IC.ARL.WECO) 1L else 2L))
      
      for (j in 1:length(shifts)) {
        runs.rules.curve[j] = get.ARL(if (IC.ARL == IC.ARL.WECO) 0L else 1L, k = NA, IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 4L, c(shifts[j], NA))
        shewhart.curve[j]   = get.ARL(2L, k = NA, IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 4L, c(shifts[j], NA))
        
        for (i in 1:length(ks)) {
          cusum.curves[i, j] = get.ARL(3L, k = ks[i], IC.ARL, steady.state.size, sigma.new = sigma.new, shift.type = 4L, c(shifts[j], NA))
        }
      }
      
      xmin = min(log(shifts))
      xmax = max(log(shifts))
      ymin = 0.0
      ymax = if (plot.capped) IC.ARL*1.2 else max(max(runs.rules.curve, shewhart.curve, cusum.curves)*1.5, IC.ARL*2.0)
      
      plot(log(shifts), runs.rules.curve, xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = bquote(log("shift")), ylab = "ARL",
           xaxt = "n",
           type = "o", lty = 1, pch = 1, col = "black", 
           main = c("Sustained shift",
                    (if (IC.ARL == IC.ARL.WECO) "WE vs. Shewhart vs. CUSUM" else "NE vs. Shewhart vs. CUSUM")))
      axis(side = 1, at = log(shifts), labels = sapply(shifts, function(s) paste("log(", s, ")", sep = "")))
      
      lines(log(shifts), shewhart.curve, type = "o", lty = 1, pch = 0, col = "red")
      sapply(1:length(ks), function(i) lines(log(shifts), as.vector(cusum.curves[i, ]), type = "o", lty = i, pch = 2, col = "blue"))
      abline(h = IC.ARL, lty = 2)
      
      legend("topright", legend = c("IC ARL", if (IC.ARL == IC.ARL.WECO) "WE" else "NE", "Shewhart",
                                    sapply(1:length(ks), function(i) paste("CUSUM(k=", ks[i], ")", sep = ""))),
             pch = c(NA, 1, 0, rep(2, length(ks))),
             lty = c(2, 1, 1, 1:length(ks)),
             col = c("black", "black", "red", rep("blue", length(ks))), bg = "white", pt.cex = 1, cex = 0.7)
    }
    
    grDevices::dev.off()
  }
}