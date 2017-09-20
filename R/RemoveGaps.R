#' @title Remove gaps from intervals 
#' @description Check whether a set of intervals (ivals) contains gaps (given as a second set of intervals). If so, either remove or split the original intervals.
#' @name RemoveGaps
#' @param gaps Dataframe generated with \code{\link{gapsts}} containing gaps in waterlevel time series
#' @param ivals Data frame of intervals that have to be corrected for gaps (typically dry times DTS or inundation times ITs)
#' @param method Method to remove gaps. "All": every interval containing (part of) a gap is removed. "Split": intervals are split into smaller intervals before and after the gap. "None" nothing is done 
#' @export RemoveGaps
#' 

RemoveGaps <- function(gaps, ivals, method=c("All", "Split")){
ivals <- ivals[,c(1,2)]
  method <- match.arg(method)
  
if (method == "All"){
  # All time intervals that contain a gap are removed from list of ITs and DTs
  for (i in (1:dim(gaps)[1])){
    gapi <- gaps[i,] 
    
    # remove all time intervals that encompass a gap
    if (any(jtest <- ivals$t1<gapi$t1&ivals$t2>gapi$t2)){
      j <- which(jtest)
      ivals <- ivals[-j,]   # remove ITs that start before beginning of gap and end after end of gap
    }

        # remove time intervals that end or start exactly at the beginning or end of a gap
    ivals <- ivals[!(ivals$t1==gapi$t1|ivals$t2==gapi$t2),]   
  }
} else if (method == "Split"){
  for (i in (1:dim(gaps)[1])){
    gapi <- gaps[i,]
    # split all time intervals that encompass a gap
    if (any(jtest <- ivals$t1<gapi$t1&ivals$t2>gapi$t2)){
      j <- which(jtest)
      ivals <-  rbind(ivals[1:j-1, ], 
                    data.frame(t1 = ivals$t1[j], t2= gapi$t1 ), # start ivals until start extended gaps
                    data.frame(t1 = gapi$t2, t2= ivals$t2[j] ), # end extended gaps until end of ivals
                    ivals[j+1:nrow(ivals),]) # if j is last row of ivals, NA row introduced
      ivals <- na.omit(ivals)
    }
    
    # remove gaps that are exactly at the beginning of a IT or DT
    if (any(jtest <- ivals$t1==gapi$t1)){
      j <- which(jtest)
      ivals$t1[j] <- gapi$t2
    }
    
    # remove gaps that are exactly at the end of a IT or DT
    if (any(jtest <- ivals$t2==gapi$t2)){
      j <- which(jtest)
      ivals$t2[j] <- gapi$t1
    }

  }
  # ivals$dt <- difftime(ivals$t1, ivals$t2, units=unit)
}  
#
return(ivals)
}
