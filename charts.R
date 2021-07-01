# (C) Michael Pokojovy and J. Marcus Jobe (2020)

## Two-Sided Shewhart Charts with Supplementary Runs Rules

# Compound rule
CRule <- function(mask) {
  ERule <- function(index = 1L) {
    switch(index,
           "1" = c(c(1, 1,  -Inf, -3.00), c(1, 1, 3.00,  Inf)), # 1
           "2" = c(c(2, 3, -3.00, -2.00), c(2, 3, 2.00, 3.00)), # 2
           "3" = c(c(4, 5, -3.00, -1.00), c(4, 5, 1.00, 3.00)), # 3
           "4" = c(c(8, 8, -3.00, -0.00), c(8, 8, 0.00, 3.00)), # 4
           "5" = c(c(2, 2, -3.00, -2.00), c(2, 2, 2.00, 3.00)), # 5
           "6" = c(c(5, 5, -3.00, -1.00), c(5, 5, 1.00, 3.00)), # 6
           "7" = c(c(1, 1,  -Inf, -3.09), c(1, 1, 3.09,  Inf)), # 7
           "8" = c(c(2, 3, -3.09, -1.96), c(2, 3, 1.96, 3.09)), # 8
           "9" = c(c(8, 8, -3.09, -0.00), c(8, 8, 0.00, 3.09))) # 9
  }
  
  res = unlist(sapply(X = mask, FUN = ERule, simplify = FALSE))

  return(matrix(res, ncol = 4, byrow = TRUE))
}

Mask2Index <- function(mask) {
  mask = as.integer(sort(mask))
  
  if (identical(mask, c(1L))) {
    1L
  } else if (identical(mask, c(7L))) {
    2L
  } else if (identical(mask, c(1L, 2L))) {
    3L
  } else if (identical(mask, c(7L, 8L))) {
    4L
  } else if (identical(mask, c(1L, 5L))) {
    5L
  } else if (identical(mask, c(1L, 3L))) {
    6L
  } else if (identical(mask, c(1L, 4L))) {
    7L
  } else if (identical(mask, c(7L, 9L))) {
    8L
  } else if (identical(mask, c(1L, 6L))) {
    9L
  } else if (identical(mask, c(1L, 2L, 3L))) {
    10L
  } else if (identical(mask, c(1L, 5L, 6L))) {
    11L
  } else if (identical(mask, c(1L, 2L, 4L))) {
    12L
  } else if (identical(mask, c(7L, 8L, 9L))) {
    13L
  } else if (identical(mask, c(1L, 3L, 4L))) {
    14L
  } else if (identical(mask, c(1L, 4L, 5L, 6L))) {
    15L
  } else if (identical(mask, c(1L, 2L, 3L, 4L))) {
    16L
  } else {
    NA
  }
}

Index2Mask <- function(index) {
  masks = list(c(1L),
               c(7L),
               c(1L, 2L),
               c(7L, 8L),
               c(1L, 5L),
               c(1L, 3L),
               c(1L, 4L),
               c(7L, 9L),
               c(1L, 6L),
               c(1L, 2L, 3L),
               c(1L, 5L, 6L),
               c(1L, 2L, 4L),
               c(7L, 8L, 9L),
               c(1L, 3L, 4L),
               c(1L, 4L, 5L, 6L),
               c(1L, 2L, 3L, 4L))
  
  return(masks[[index]])
}

Mask2IC.ARL <- function(mask) {
  index = Mask2Index(mask)
  
  ARLs = c(370.40, 499.62, 225.44, 239.75, 278.03, 166.05, 152.73, 170.41, 349.38, 132.89, 266.82, 122.05, 126.17, 105.78, 133.21, 91.75)
  
  return(ARLs[index])
  
}

check_CompoundChart <- function(object) {
  errors <- character()

  if (length(object@mask) == 0) {
    msg <- paste("mask should be non-empty")
    errors <- c(errors, msg)
  } else {
    if (max(object@mask) > 9L) {
      msg <- paste("Entries of mask should be 1...9")
      errors <- c(errors, msg)
    }
  }

  if (length(object@mask) > 9) {
    msg <- paste("Length of mask should not exceed 9")
    errors <- c(errors, msg)
  }

  if (length(errors) == 0) TRUE else errors
}

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

CreateCompoundChart <- function(mask = c(1L)) {
  object = new("CompoundChart", mask = mask)
  
  object@time    = 0L
  object@stopped = FALSE
  object@rule    = CRule(mask)
  object@state   = sapply(1:nrow(object@rule), function(i) logical(length = object@rule[i, 2]), simplify = FALSE) 

  return(object)
}

UpdateCompoundChart <- function(object, x) {
  if (object@stopped) {
    return(object)
  } else {
    object@time = object@time + 1L
    
    for (i in 1:nrow(object@rule)) {
      flag = (object@rule[i, 3] < x) && (x < object@rule[i, 4])
      len  = length(object@state[i]) 
      
      object@state[[i]] = c(object@state[[i]][-1], flag)
      
      if (sum(object@state[[i]]) == object@rule[i, 1]) {
        object@stopped = TRUE
      }
    }
    
    return(object)
  }
}

## Two-Sides Shewhart Chart

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

CreateShewhartChart <- function(IC.ARL = 370, CL = qnorm(0.5/IC.ARL, lower.tail = FALSE)) {
  object = new("ShewhartChart", IC.ARL = IC.ARL)
  
  object@CL      = if (is.na(CL)) qnorm(0.5/IC.ARL, lower.tail = FALSE) else CL
  object@time    = 0L
  object@stopped = FALSE
  
  return(object)
}

UpdateShewhartChart <- function(object, x) {
  if (object@stopped) {
    return(object)
  } else {
    object@time    = object@time + 1L
    object@stopped = (abs(x) > object@CL)
    
    return(object)
  }
}

## Two-Sided CUSUM Chart

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

CreateCUSUMChart <- function(IC.ARL = 370, k = 0.5, CL = spc::xcusum.crit(k, L0 = IC.ARL, sided = "two")) {
  object = new("CUSUMChart", IC.ARL = IC.ARL, CL = CL)
  
  object@CL      = if (is.na(CL)) spc::xcusum.crit(k, L0 = IC.ARL, sided = "two") else CL
  object@k       = k
  object@time    = 0L
  object@stopped = FALSE
  object@state   = c(0.0, 0.0)
  
  return(object)
}

UpdateCUSUMChart <- function(object, x) {
  if (object@stopped) {
    return(object)
  } else {
    object@time  = object@time + 1L
    object@state = c(max(0.0, object@state[1] + (x - object@k)), min(0.0, object@state[2] + (x + object@k)))
    
    object@stopped = (object@state[1] > object@CL) || (object@state[2] < -object@CL)
    
    return(object)
  }
}

## Nelson's chart

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

UpdateNelsonChart <- function(object, x) {
  if (object@stopped) {
    return(object)
  } else {
    object@time  = object@time + 1L
    
    object@cnt.rule1  = as.integer(abs(x) > 3.0)
    object@cnt.rule2a = if (x <= 0.0) object@cnt.rule2a + 1L else 0L
    object@cnt.rule2b = if (x >= 0.0) object@cnt.rule2b + 1L else 0L
    
    if (!is.na(object@x.prev[3])) {
      object@cnt.rule3a = if ((x - object@x.prev[3]) <= 0) object@cnt.rule3a + 1L else 0L
      object@cnt.rule3b = if ((x - object@x.prev[3]) >= 0) object@cnt.rule3b + 1L else 0L
    }
    
    x   = c(object@x.prev, x)
    dx2 = diff(x, lag = 1L, differences = 2L)
    if (sum(is.na(dx2)) == 0L) {
      object@cnt.rule4 = if (prod(dx2) <= 0.0) {if (object@cnt.rule4 > 0L) object@cnt.rule4 + 1L else 4L} else 0L
    } 
    object@x.prev = x[2:4]
    x = x[4]
      
    object@cnt.rule7 = if (abs(x) <= 1.0) object@cnt.rule7 + 1L else 0L
    object@cnt.rule8 = if (abs(x)  > 1.0) object@cnt.rule8 + 1L else 0L
    
    object@state.rule5a = c(object@state.rule5a[2:3], (x < -2.0))
    object@state.rule5b = c(object@state.rule5b[2:3], (x >  2.0))
    object@state.rule6a = c(object@state.rule6a[2:5], (x < -1.0))
    object@state.rule6b = c(object@state.rule6b[2:5], (x >  1.0))
    
    object@stopped = object@stopped || (object@cnt.rule1  >= 1L)
    object@stopped = object@stopped || (object@cnt.rule2a >= 9L)
    object@stopped = object@stopped || (object@cnt.rule2b >= 9L)
    object@stopped = object@stopped || (object@cnt.rule3a >= 6L)
    object@stopped = object@stopped || (object@cnt.rule3b >= 6L)
    object@stopped = object@stopped || (object@cnt.rule4  >= 14L)
    object@stopped = object@stopped || (object@cnt.rule7  >= 15L)
    object@stopped = object@stopped || (object@cnt.rule8  >= 8L)
    
    object@stopped = object@stopped || (sum(object@state.rule5a) == 2L)
    object@stopped = object@stopped || (sum(object@state.rule5b) == 2L)
    object@stopped = object@stopped || (sum(object@state.rule6a) == 4L)
    object@stopped = object@stopped || (sum(object@state.rule6b) == 4L)
    
    return(object)
  }
}