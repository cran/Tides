gapsts <- function(ts,            	# array of times, consisting of different continuous subseries seperated by large gaps
                   dtMax,        	# maximum time interval in a continuous series, or equivalently minimum interval to be characterized as gap.
                   unit = "mins",  	# unit of dtMax; used when ts belongs to class POSIXt
                   shiftbegin = FALSE	# add time step between fore-last and last measurement before gap to start time (t0) of gap	
)
{
  if (!inherits(ts, "POSIXt")) {
    timediffs <- ts[1:(length(ts) - 1)] - ts[2:(length(ts))]
  }
  else {
    timediffs <- difftime(ts[1:(length(ts) - 1)], ts[2:(length(ts))], 
                          units = unit)
  }
  
  if (!any(timediffs < - dtMax)) return(NULL)
  #Select gaps > dtMax in a timeseries ts
  gaps <- ts[c(timediffs < -dtMax,FALSE)]
  gaps <- data.frame(t1 = gaps)
  gaps$t2 <- ts[match(gaps$t1,ts)  + 1]
  gaps$n <- 1:dim(gaps)[1]
  
  
  if (shiftbegin){
    # if (length(ts) < 3) return("Time series consists of less than 3 data points")
    
    t1_i <- match(gaps$t1, ts)
    if (t1_i[1] == 1){
      t1_i[1] <- which(Mod(timediffs) <= dtMax)[1] + 1 # ADJUSTED! # 3
      # if first data point is beginning of gap, dt cannot be estimated from the former continuous series, but is estimated from the next continuous series
      # ADJUSTED! it's also possible that the second difftime > dtMax!
      # take the first difftime < dtMax (+1, since in the next step one is subtracted)
      warning("First data point is beginning of gap. To shift t0, dt is estimated from next continuous series")
    }
    gaps$t1 <- gaps$t1 + Mod(timediffs[t1_i-1])
    
  }
  
  # CHANGED!
  # moved this part from above Shiftbegin, otherwise dt is not calculated with the adjusted gaps$t1 ???
  if (!inherits(ts,"POSIXt")){
    gaps$dt <- gaps$t2 - gaps$t1
  } else {
    gaps$dt <- difftime(gaps$t2,gaps$t1,units=unit)
  }
  
  return(gaps) #Data frame with the initial time, end time and time difference (unit = unit) of each interval > dtMax
}



