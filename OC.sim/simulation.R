# (C) Michael Pokojovy and J. Marcus Jobe (2020)

do.simulation <- function(chart.index = 0L, k = NA, IC.ARL = NA, steady.state.size = 0L, nrep = 1E4L, nmax = 1E6L, 
                          sigma.new = 1.0, shift.type = 0L, par = NULL) {
  if (shift.type == 0L) {
    # sustained: par[1]
    shift.fun = function(n) 0.0
    par = c(NA, NA)
  } else if (shift.type == 1L) {
    # slope: par[1]
    shift.fun = function(n) par[1]*n
    par = c(par[1], NA)
  } else if (shift.type == 2L) {
    # amplitude: par[1], period: par[2]
    shift.fun = function(n) par[1]*sin(2*pi*n/par[2])
    par = par[1:2]
  } else if (shift.type == 3L) {
    # amplitude: par[1], period: par[2]
    shift.fun = function(n) ifelse((n - 1) %% par[2] < par[2]/2, par[1], -par[1])
    par = par[1:2]
  } else if (shift.type == 4L) {
    # shift: par[1]
    shift.fun = function(n) par[1]
    par = c(par[1], NA)
  } else {
    par = c(NA, NA)
  }

  if (chart.index == 2L) {
    CL = qnorm(0.5/IC.ARL, lower.tail = FALSE)
  } else if (chart.index == 3L) {
    CL = spc::xcusum.crit(k, L0 = IC.ARL, sided = "two")
  }
  
  ARL  = 0.0
  ARL2 = 0.0
  
  surv.rate    = 0.0
  overrun.rate = 0.0
  
  iteration <- function(chart.index, k = NA, CL = NA, IC.ARL = NA, steady.state.size = steady.state.size, 
                        sigma.new = sigma.new, shift.fun = shift.fun) {
    if (chart.index == 0L) {
      setClass("CompoundChart",
               slots = c(
                 mask    = "integer",
                 time    = "integer",
                 stopped = "logical",
                 rule    = "matrix",
                 state   = "list"
               ),
               prototype = list(
                 mask    = integer(),
                 time    = 0L,
                 stopped = FALSE,
                 rule    = matrix(),
                 state   = list()
               ),
               validity = check_CompoundChart
      )
      chart = CreateCompoundChart(c(1L, 2L, 3L, 4L))
    } else if (chart.index == 1L) {
      setClass("NelsonChart",
               slots = c(
                 time       = "integer",
                 stopped    = "logical",
                 cnt.rule1  = "integer",
                 cnt.rule2a = "integer",
                 cnt.rule2b = "integer",
                 cnt.rule3a = "integer",
                 cnt.rule3b = "integer",
                 cnt.rule4  = "integer",
                 cnt.rule7  = "integer",
                 cnt.rule8  = "integer",
                 state.rule5a = "logical",
                 state.rule5b = "logical",
                 state.rule6a = "logical",
                 state.rule6b = "logical",
                 x.prev = "numeric"
               ),
               prototype = list(
                 time       = 0L,
                 stopped    = FALSE,
                 cnt.rule1  = 0L,
                 cnt.rule2a = 0L,
                 cnt.rule2b = 0L,
                 cnt.rule3a = 0L,
                 cnt.rule3b = 0L,
                 cnt.rule4  = 0L,
                 cnt.rule7  = 0L,
                 cnt.rule8  = 0L,
                 state.rule5a = rep(FALSE, 3),
                 state.rule5b = rep(FALSE, 3),
                 state.rule6a = rep(FALSE, 5),
                 state.rule6b = rep(FALSE, 5),
                 x.prev = as.numeric(rep(NA, 3L))
               )
      )
      chart = new("NelsonChart")
    } else if (chart.index == 2L) {
      setClass("ShewhartChart",
               slots = c(
                 IC.ARL  = "numeric",
                 CL      = "numeric",
                 time    = "integer",
                 stopped = "logical"
               ),
               prototype = list(
                 IC.ARL  = integer(),
                 CL      = integer(),
                 time    = 0L,
                 stopped = FALSE
               )
      )
      chart = CreateShewhartChart(IC.ARL = IC.ARL, CL = CL)
    } else {
      setClass("CUSUMChart",
               slots = c(
                 IC.ARL  = "numeric",
                 k       = "numeric",
                 CL      = "numeric",
                 time    = "integer",
                 stopped = "logical",
                 state   = "numeric"
               ),
               prototype = list(
                 IC.ARL  = integer(),
                 k       = numeric(),
                 CL      = integer(),
                 time    = 0L,
                 stopped = FALSE,
                 state   = c(0.0, 0.0)
               )
      )
      chart = CreateCUSUMChart(IC.ARL = IC.ARL, k = k, CL = CL)
    }

    survived = TRUE

    if (steady.state.size > 0) {
      for (RL in 1:steady.state.size) {
        if (chart.index == 0) {
          chart = UpdateCompoundChart(chart, rnorm(1))
        } else if (chart.index == 1) {
          chart = UpdateNelsonChart(chart, rnorm(1))
        } else if (chart.index == 2) {
          chart = UpdateShewhartChart(chart, rnorm(1))
        } else {
          chart = UpdateCUSUMChart(chart, rnorm(1))
        }

        if (chart@stopped) break
      }
    }

    RL = 0L

    if (chart@stopped) {
      survived = FALSE
    } else {
      for (i in 1L:nmax) {
        RL = RL + 1L

        if (chart.index == 0) {
          chart = UpdateCompoundChart(chart, sigma.new*rnorm(1) + shift.fun(RL))
        } else if (chart.index == 1) {
          chart = UpdateNelsonChart(chart, sigma.new*rnorm(1) + shift.fun(RL))
        } else if (chart.index == 2) {
          chart = UpdateShewhartChart(chart, sigma.new*rnorm(1) + shift.fun(RL))
        } else {
          chart = UpdateCUSUMChart(chart, sigma.new*rnorm(1) + shift.fun(RL))
        }

        if (chart@stopped) break
      }
    }
    
    overrun = (RL >= nmax)

    return(c(survived, overrun, RL, RL^2, RL*(!overrun), RL^2*(!overrun))/nrep)
  }
  
  if (parallel.flag) {
    if (chart.index == 0L) {
      func.list = c("CRule", "check_CompoundChart", "CreateCompoundChart", "UpdateCompoundChart")
    } else if (chart.index == 1L) {
      func.list = c("UpdateNelsonChart")
    } else if (chart.index == 2L) {
      func.list = c("CreateShewhartChart", "UpdateShewhartChart")
    } else {
      func.list = c("CreateCUSUMChart", "UpdateCUSUMChart")
    }

    sim.res = foreach(ind = 1:nrep, .inorder = FALSE, .export = func.list, .combine = "+") %dopar% 
                      iteration(chart.index = chart.index, k = k, CL = CL, IC.ARL = IC.ARL, steady.state.size = steady.state.size, 
                                sigma.new = sigma.new, shift.fun = shift.fun)
  } else {
    sim.res = foreach(ind = 1:nrep, .inorder = TRUE, .combine = "+") %do% 
                      iteration(chart.index = chart.index, k = k, CL = CL, IC.ARL = IC.ARL, steady.state.size = steady.state.size, 
                                sigma.new = sigma.new, shift.fun = shift.fun)
  }
  
  surv.rate    = sim.res[1]
  overrun.rate = sim.res[2]
  ARL          = sim.res[3] # ARL over ALL streams
  ARL2         = sim.res[4] # Avg. RL^2 over ALL streams
  ARL0         = sim.res[5] # ARL over non-truncated streams only
  ARL20        = sim.res[6] # Avg. RL^2 over non-truncated streams only
  
  overrun.rate = overrun.rate/surv.rate
  
  nrep.cor  = nrep*surv.rate
  nrep.cor0 = nrep*surv.rate*(1 - overrun.rate)
  
  ARL0  = if (nrep.cor0 <= 0) 0.0 else ARL0*(nrep/nrep.cor0)
  ARL20 = if (nrep.cor0 <= 1) 0.0 else ARL20*(nrep/(nrep.cor0 - 1))
  
  ARL  = ARL*(nrep/nrep.cor)
  ARL2 = ARL2*(nrep/(nrep.cor - 1))
  
  VRL0    = max(0.0, (ARL20 - (nrep.cor0/(nrep.cor0 - 1))*ARL0^2))
  std.RL0 = sqrt(max(0.0, VRL0))
  
  VRL    = max(0.0, (ARL2 - (nrep.cor/(nrep.cor - 1))*ARL^2))
  std.RL = sqrt(max(0.0, VRL))
  
  res = matrix(c(chart.index, k, IC.ARL, steady.state.size, ARL, VRL, std.RL, ARL0, VRL0, std.RL0,
                 surv.rate, overrun.rate, sigma.new, shift.type, par), nrow = 1L)
  colnames(res) = c("chart.index", "k", "IC.ARL", "steady.state.size", "ARL", "VRL", "std.RL", "ARL0", "VRL0", "std.RL0",
                    "surv.rate", "overrun.rate", "sigma.new", "shift.type", "par.1", "par.2")
  return(res)
}