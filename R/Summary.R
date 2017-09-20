
#' @title Summary of tidal characteristics
#' 
#' @description Outputs list of summary statistics of a Tides oabject
#' 
# #' @usage S3 method for class 'Tides'
#' 
#' @param object Tides object (e.g. the result of \code{\link{TidalCharacteristics}})
#' @param ... Not used (Added for S3 method compatibility)
#' @author {Lennert Schepers <Lennert.Schepers@uantwerp.be>, Tom Cox <tom.cox@uantwerp.be>}
#' @export summary.Tides
#' 
#' @return 
#' A list containing: 
#' \itemize{
#' \item{IFsum:}{ Inundation frequency: see IF(). The proportion of high water tides higher than h0. A warning will be displayed when the inundation frequency is 100\%.}
#' \item{nIndsum:}{ Inundations during time span, when the waterlevel > h0.}
#' \item{IHsum:}{ Average inundation height. A warning will be displayed when the inundation frequency is 100\%.}
#' \item{IHcsum:}{ Average inundation height (per cycle).}
#'
#' \item{Tunitssum:}{ time units}
#' \item{ITsum:}{ Average inundation time, in Tunitssum A warning will be displayed when h0 is never inundated}
#' \item{ITCsum:}{ Average inundation time (per cycle), in Tunitssum}
#' \item{ITMsum:}{ Maximal inundation time, in Tunitssum}
#' \item{DTsum:}{ Average dry time, in Tunitssum}
#' \item{DTCsum:}{ Average dry time (per cycle), in Tunitssum}
#' \item{DTMsum:}{ Maximal dry time, in Tunitssum}
#'
#' \item{MHWsum:}{ Average high water. Note that the calculated HW and LW are always >= h0!}
#' \item{MLWsum:}{ Average low water. Note that the calculated HW and LW are always >= h0!}
#' \item{TRsum:}{ Average tidal range. Note that the calculated HW and LW are always >= h0!}
#'
#' \item{nTCsum:}{ number of (tidal) cycles}
#' \item{nTCFsum:}{ number of full tidal cycles (used to measure averages per cycle)}
#'
#' \item{nGsum:}{ number of gaps}
#' \item{GTsum:}{ total gaps time in mins}
#' \item{nTSsum:}{ number of continuous timeseries}
#' \item{TTNoGapsum:}{ total continuous timeseries time (without gaps, also tidal phases before and after gap are not included)}
#'
#' \item{ITTsum:}{ total inundation time (without gaps, and also tidal phases before and after gap are not included)}
#' \item{DTTsum:}{ total dry time (without gaps, and also tidal phases before and after gap are not included)}

#' \item{IPsum:}{ proportion of total time inundated (without gaps, and also tidal phases before and after gap are not included)}
#' \item{DPsum:}{ DTTsum/TSsum # proportion of total time dry (without gaps, and also tidal phases before and after gap are not included)}
#' }

summary.Tides<- function(object, ...){
  x <- object   # internal alias
  nTCsum <- x$Ncycles # number of tidal cycles # updated method: see TidalCharacteristics()
  nTCFsum <- x$Nfullcycles
  if (is.null(x$DTs$dt)) {
    warning("Inundation frequency: 100% \n")
    IFsum <- 100
    nIndsum <- 0
  } else {
    IFsum <- x$IF
    nIndsum <- x$IF*x$Ncycles/100
  }
  
  if (is.null(x$DTs$dt)) {
    IHsum <- mean(x$HL$h-x$HL$h0)
    IHcsum <- sum(x$HL$h[x$HL$HL=="H"]-x$HL$h0[x$HL$HL=="H"])/x$Ncycles
    warning("WARNING not reliable (IF = 100%). Average inundation height: ", IHsum ,"\n")
  } else {
    IHsum <- mean(x$HL$h[x$HL$HL=="H"]-x$HL$h0[x$HL$HL=="H"])
    IHcsum <- sum(x$HL$h[x$HL$HL=="H"]-x$HL$h0[x$HL$HL=="H"])/x$Nfullcycles
  }
  
  # time units  
  Tunitssum <- x$Tunit
  
  if (is.null(x$ITs$dt)) {
    warning("Average inundation time: Site never inundated")
    ITsum <- 0
    ITCsum <- 0
    ITMsum <- 0
  } else {
    ITsum <- mean(x$ITs$dt)
    ITCsum <- sum(x$ITs$dt)/x$Ncycles
    ITMsum <- max(x$ITs$dt)
  }
  
  if (is.null(x$DTs$dt)) {
    warning("Average dry time: Site never falls dry")
    DTsum <- 0
    DTCsum <- 0
    DTMsum <- 0
  } else {  
    DTsum <- mean(x$DTs$dt)
    DTCsum <- sum(x$DTs$dt)/x$Ncycles
    DTMsum <- max(x$DTs$dt)
  }
  
  MHWsum <- mean(x$HL$h[x$HL$HL=="H"]) # Note that the calculated high water is always >= h0!
  MLWsum <- mean(x$HL$h[x$HL$HL=="L"]) # Note that the calculated low water is always >= h0!
  
  # Is it useful to calculate the mean tidal range since all values <h0 are replaced by h0?
  # yes! If you put h0 lower than lowest value, you can calculate the 'real' mean tidal range in a normal tidal cycle
  TRsum <- MHWsum-MLWsum # mean tidal range. Note that the calculated HW and LW are always >= h0!
  
  
  
  if (is.null(x$gaps)) {
    warning("There were no gaps in the time series")
    GTsum <- 0
    nGsum <- 0
  } else {
    GTsum <- sum(unclass(x$gaps$dt)) # in mins
    # # TO DO: add removed gaptime? (for drytime and inundation time??)
    #  x$h$time[is.element(x$h$c(x$gaps$N-1,x$gaps$N,x$gaps$N+1)))
    # GTDelsum <- 
    nGsum <- max(x$gaps$n)
  }
  
  nTSsum <- max(x$h$n) # number of continuous timeseries
  
  
  ITTsum <- sum(x$ITs$dt) # total inundation time (without gaps, and also tidal phases before and after gap are not included)
  DTTsum <- sum(x$DTs$dt) # total dry time (without gaps, and also tidal phases before and after gap are not included)
  TTNoGapsum <- ITTsum + DTTsum # total continuous timeseries time (without gaps, and also tidal phases before and after gap are not included)
  
  IPsum <-  as.double(ITTsum, units='secs')/as.double(TTNoGapsum, units='secs') # proportion of total time inundated (without gaps, and also tidal phases before and after gap are not included)
  DPsum <- as.double(DTTsum, units='secs')/as.double(TTNoGapsum, units='secs') # proportion of total time dry (without gaps, and also tidal phases before and after gap are not included)
  
  
  return(list(IFsum = IFsum, nIndsum= nIndsum, IHsum=IHsum, IHcsum=IHcsum, 
              Tunitssum=Tunitssum, ITsum=ITsum, ITCsum=ITCsum, ITMsum=ITMsum,
              DTsum=DTsum, DTCsum=DTCsum, DTMsum=DTMsum,
              MHWsum=MHWsum, MLWsum=MLWsum, TRsum=TRsum, nTCsum=nTCsum,
              nTCFsum=nTCFsum,
              nGsum=nGsum, GTsum=GTsum, nTSsum=nTSsum, TTNoGapsum=TTNoGapsum,
              ITTsum=ITTsum, DTTsum=DTTsum, IPsum=IPsum, DPsum=DPsum))
  
}

