# (C) Michael Pokojovy and J. Marcus Jobe (2020)

do.simulation <- function(chart.index = 0L, k = NA, IC.ARL = NA, nrep = 1E4L, nmax = 1E6L, new.sigma = 1.0, shift.type = 0L, par = NULL) {
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
    shift.fun = if (period > 3) function(n) par[1]*sin(2*pi*n/par[2]) else function(n) par[1]*sin(2*pi*n/par[2])/(sqrt(3)/2)
    par = par[1:2]
  } else if (shift.type == 3L) {
    # amplitude: par[1], period: par[2]
    shift.fun = function(n) ifelse((n - 1) %% par[2] < par[2]/2, par[1], -par[1])
    par = par[1:2]
  } else {
    par = c(NA, NA)
  }

  if (chart.index == 2L) {
    CL = qnorm(0.5/IC.ARL, lower.tail = FALSE)
  } else if (chart.index == 3L) {
    CL = spc::xcusum.crit(k, L0 = IC.ARL, sided = "two")
  }
  
  iteration <- function(chart.index, k = NA, CL = NA, IC.ARL = NA, new.sigma = new.sigma, shift.fun = shift.fun) {
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

    RL = 0L

    for (i in 1:nmax) {
      RL = RL + 1L
      
      if (chart.index == 0) {
        chart = UpdateCompoundChart(chart, new.sigma*rnorm(1) + shift.fun(RL))
      } else if (chart.index == 1) {
        chart = UpdateNelsonChart(chart, new.sigma*rnorm(1) + shift.fun(RL))
      } else if (chart.index == 2) {
        chart = UpdateShewhartChart(chart, new.sigma*rnorm(1) + shift.fun(RL))
      } else {
        chart = UpdateCUSUMChart(chart, new.sigma*rnorm(1) + shift.fun(RL))
      }
      
      if (chart@stopped) break
    }
    
    return(RL)
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

    sim.res = foreach(ind = 1:nrep, .inorder = FALSE, .export = func.list, .combine = "c") %dopar% 
                      iteration(chart.index = chart.index, k = k, CL = CL, IC.ARL = IC.ARL, new.sigma = new.sigma, shift.fun = shift.fun)
  } else {
    sim.res = foreach(ind = 1:nrep, .inorder = TRUE, .combine = "c") %do% 
                      iteration(chart.index = chart.index, k = k, CL = CL, IC.ARL = IC.ARL, new.sigma = new.sigma, shift.fun = shift.fun)
  }
  
  return(sim.res)
}