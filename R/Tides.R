
#' @name TidalCharacteristics
#' @aliases TidalCharacteristics
#' @title Calculate tidal characteristics
#' @description Calculates the characteristics of observed tidal water levels. Wrapper of the functions \code{\link{extrema}}, \code{\link{IT}} and \code{\link{IF}}. Also works on time series with gaps.
#' @param h Water level time series. data frame with time and h column
#' @param h0 Reference level, either single valued or vector with dimension corresponding to h
#' @param T2 'Lower' bound on half the quasi period, but higher than expected stagnant phase; default = 5h
#' @param hoffset Offset level, to prevent spurious maxima generation due to small fluctuations
#' @param filtconst Filtering constant for smoothing the time series
#' @param dtMax Maximum accepted time interval in a continuous series. Bigger time intervals are considered to be gaps
#' @param unit Unit of dtMax, Tavg
#' @param Tavg Average period of time series
#' @param removegaps Method to remove gaps in time series from inundation times and dry times. See \code{\link{RemoveGaps}}
#' @return An object of class \code{Tides}, i.e. a list containing:
#' \itemize{
#'  \item{HL}{Data frame with extrema}
#'  \item{h }{original water level data frame with additional attributes}
#'  \item{gaps}{a data frame containing start and end times of gaps in the series}
#'  \item{IF}{inundation frequency of the reference level}
#'  \item{ITs}{inundation times at the reference level}
#'  \item{DTs}{dry times at the reference level}
#'  \item{h0}{reference level}
#'  \item{N}{Total number of cycles in time span}
#' }
#' @seealso \code{\link{extrema}}, \code{\link{IT}}, \code{\link{plot.Tides}}
#' @author Tom Cox <tom.cox@uantwerp.be>, Lennert Schepers <lennert.schepers@uantwerp.be>
#' @keywords utilities
#' @examples 
#' TC <- TidalCharacteristics(waterlevels, filtconst=10,hoffset=1)
#'  TC
#'  plot(TC)
#'  summary(TC)




TidalCharacteristics <- function (	h,  		#(Water level) time series. data frame with time and h column
					h0 = h$h0, 	#Reference level, either single valued or vector with dimension corresponding to h
					T2 = 5*60*60, 	#'Lower' bound on half the quasi period, but higher than expected stagnant phase; default = 5h
					hoffset = 0,   #[Maybe obsolete]Margin on reference level, to cope with small fluctuations in the Water level time series
					filtconst = 1,	#Filtering constant for smoothing the time series
					dtMax = 15,	#maximum accepted time interval in a continuous series. Bigger time intervals are considered to be gaps
					unit = "mins",  # unit of dtMax, Tavg
					Tavg = 12.4*60, #Average period of tidal cycle 
          removegaps = c("All", "Split", "None"))
{
removegaps <- match.arg(removegaps)
if (any(!is.element("time",names(h)))) stop("h must be a data frame with columns time (POSIXt) and h (numeric)")
if (!is.element("POSIXt",class(h$time[1]))) stop("h must be a data frame with columns time (POSIXt) and h (numeric)")
if (is.null(h0)) stop("provide reference level h0")

#Calculate high and low water levels (extrema) in water level time series
output <- extrema(h=h,h0=h0,T2 = T2,filtconst=filtconst, hoffset=hoffset)
HL <- output$HL
tij <- output$h

###
#determine gaps in 'continuous' data
gaps <- gapsts(tij$time,dtMax=dtMax, shiftbegin=TRUE)
if (!is.null(gaps)) 	{gaps$N <- tij$N[match(gaps$t1,tij$time)]	# gaps$N is tidal phase at which gap starts 
			tij$n <- findInterval(tij$time,gaps$t2)+1	#n counts the continuous series of cycles.	
} else tij$n <- 1

###
# Calculate inundation times and dry times
# The functions IT and DT will determine the parts of the time series that the waterlevel is above or below the reference level
# It does not take care of potential gaps in the series. That will be done in the next step

ITDT <- IT(tij[c("time","h")],h0=tij$h0,dtMax=dtMax, hoffset = hoffset)
ITs <- ITDT$IT
if (!is.null(ITs)) ITs$N <- tij$N[match(ITs$t1,tij$time)] # ITs$N = tidal phase at which inundation starts
  
DTs <- ITDT$DT
if (!is.null(DTs)) DTs$N <- tij$N[match(DTs$t1,tij$time)] # DTs$N = tidal phase at which dry time starts


#
# Remove gaps from ITs and DTs
#

#
# OLD CODE
#

# When gaps are present in the time series, remove all broken cycles +- 1 cycle number to be sure, i.e. in which a gap of data exists
# if (!is.null(gaps)) {
#  ITs <- ITs[!is.element(ITs$N, c(gaps$N-1,gaps$N,gaps$N+1)),]
#  DTs <- DTs[!is.element(DTs$N,gaps$N)&!is.element(DTs$N,gaps$N+1)&!is.element(DTs$N,gaps$N - 1),]
#} 

if (!is.null(gaps)) {
  ITs_nogaps <- RemoveGaps(gaps, ivals=ITs, method=removegaps)
  ITs <- cbind(ITs_nogaps, ITs[match(ITs_nogaps$t1, ITs$t1), c("n", "dt", "N")])
  DTs_nogaps <- RemoveGaps(gaps, ivals=DTs, method=removegaps)
  DTs <- cbind(DTs_nogaps, DTs[match(DTs_nogaps$t1, DTs$t1), c("n", "dt", "N")])
  
} 



# ADDED Lennert Schepers:
#---------------
# the code above assumes that 1 inundation/dry period lasts only  tide!
# solution: delete time of gaps +-1N from inundation/dry period

# if (!is.null(gaps)) { 
# if(!is.null(DTs)) DTs <- CutExtendedGaps(tij=tij, gaps=gaps, ivals=DTs)
# if(!is.null(ITs)) ITs <- CutExtendedGaps(tij=tij, gaps=gaps, ivals=ITs)
# }


#--------------------
# ADDED
# number of tidal phases that are removed for the gaps:
# unique() to remove doubles by possible successive gaps
# the first and last measurement are mostly(?) part of a separate tidal phase, so delete these 2
# we check if these are part of the deleted gaps phases: 1 and last are included in unique())
if(is.null(gaps)) Ndel <- 2 else Ndel <- length(unique(c(gaps$N-1,gaps$N,gaps$N+1, 1, max(tij$N))))


###
#Calculate total inundation frequence
if (is.null(gaps)) gapstime <- 0 else gapstime <- sum(unclass(gaps$dt))

t_total_mins <- difftime(max(tij$time,na.rm=T),min(tij$time,na.rm=T),units="mins")
t_observations_mins <- unclass(t_total_mins - gapstime)

N_observations <- floor(t_observations_mins/(Tavg))[1]


# Remark! Ncycles is calculated by dividing No gaps time/ Tavg(=input in function!) better to calculate the number of cycles?
# total tidal cycles = (total phases - removed phases for gaps -2) / 2
# two phases: H, L -> 1 tidal cycle (H is the full high water phase, L is the full low water phase, not only from HW to LW!)
# n phases -> n/2 tidal cycles
# delete first and last (incomplete) tidal stage (first and last point are separate cycle)
N_inundations <- (max(tij$N) - 2)/ 2
N_full_inundations <- (max(tij$N) - Ndel)/ 2 #counted by the number of full phases that are not deleted for gaps

IF <- IF(HL[HL$HL=="H",]$h,HL[HL$HL=="H",]$h0,N=N_observations)


TideChars <- list(HL=HL,h=tij,gaps=gaps,IF=IF,ITs=ITs,DTs=DTs,h0 = h0, t_total_mins = t_total_mins, t_observations_mins = t_observations_mins, N_observations = N_observations, N_inundations=N_inundations,N_full_inundations = N_full_inundations, Tunit=unit)
class(TideChars) <- "Tides" 
return(TideChars)
}

#
# END TidalCharacteristics
#


print.Tides <- function(x,...){
  if (is.null(x$DTs$dt)) {cat("Inundation frequency: 100% \n") 
  } else {
    cat("Inundation frequency: ", x$IF, "\n")
    cat("Inundations during time span: ", x$N_inundations, "\n")}
  
  if (is.null(x$DTs$dt)) {cat ("WARNING not reliable (IF = 100%). Average inundation height: ", mean(x$HL$h-x$HL$h0),"\n") 
  } else {
    cat("Average inundation height: ", mean(x$HL$h[x$HL$HL=="H"]-x$HL$h0[x$HL$HL=="H"]),"\n")
    cat("Average inundation height (per cycle): ", sum(x$HL$h[x$HL$HL=="H"]-x$HL$h0[x$HL$HL=="H"])/x$N_inundations,"\n")}
  
  if (is.null(x$ITs$dt)) {cat ("Average inundation time: Site never inundated \n") 
  } else {
    cat("Average inundation time: ", mean(x$ITs$dt), x$Tunit, "\n") 
    cat("Average inundation time (per cycle): ", sum(x$ITs$dt)/x$N_full_inundations, x$Tunit, "\n") # changed to Nfullcycles
    cat("Maximal inundation time: ", max(x$ITs$dt), x$Tunit, "\n")
  }
  
  
  if (is.null(x$DTs$dt)) {cat ("Average dry time: Site never falls dry \n") 
  } else {  
    cat("Average dry time: ", mean(x$DTs$dt), x$Tunit, "\n")
    cat("Average dry time (per cycle): ", sum(x$DTs$dt)/x$N_full_inundations, x$Tunit, "\n") # changed to Nfullcycles
    cat("Maximal dry time: ", max(x$DTs$dt), x$Tunit, "\n")}
  
  cat("Average high water: ", mean(x$HL$h[x$HL$HL=="H"]),"\n")
  cat("Average low water: ", mean(x$HL$h[x$HL$HL=="L"]),"\n")
  
  # ADJUSTED: Ncycles is based on arbitrary Tavg (=input data??)
  # I changed the calculation of Ncycles, see TidalCharacteristics()
  cat("Time span: ", x$t_total_mins/60/24, " days","\n")                
  cat("Observations in time span: ", x$t_observations_mins/60/24, " days or ", x$N_observations, " average tidal cycles" ,"\n")
  cat("Inundations in time span: ", floor(x$N_inundations), "tidal cycles", "\n")
  cat("Fully observed inundated tidal cycles: ", floor(x$N_full_inundations), " full tidal cycles","\n")        # ADDED
  
  if (is.null(x$gaps)) {cat("There were no gaps in the time series","\n") 
  } else {
    # ADJUSTED! number of continuous series is not x$gaps$n(=number of gaps),
    # but x$h$n
    cat("The time series consists of",max(x$h$n),"continuous sub-series","\n")}
  
}

plot.Tides <- function(x,...){
  plot(x$h$time,x$h$h,type="l",xlab="", ylab="waterlevel",...)
  lines(x$h$time,x$h$h0)
  points(x$HL$time,x$HL$h,col="red",pch=20)
  points(x$HL$time[x$HL$HL=="L"],x$HL$h[x$HL$HL=="L"],col="blue",pch=20)
}


extrema <- function(h,            #(Water level) time series. Data frame with time and h column
                    h0,           #Reference level, either single valued or vector with dimension corresponding to h
                    T2 = 5*60*60, #'Lower' bound on half the quasi period, but higher than expected stagnant phase; default = 5h
                    hoffset = 0, #Offset level, to prevent spurious maxima generation due to small fluctuations
                    filtconst = 1 #Filtering constant for smoothing the time series
)
{

h$h0 <- h0

#set all levels < h0 equal to h0
h$ho <- h$h 	#first save original waterlevels
h$h[(h$h)<=h$h0] <- h$h0[(h$h)<=h$h0]

#filter tij data, running average of filtconst (default = 1, no filtering) succesive datapoints,to
#remove small fluctuations

h$hfilt <- filter(h$h,rep(1/filtconst,filtconst))

#remove missing values due to filtering
h <- h[!is.na(h$hfilt),]

#Useful matrix to swith between [H,L] or [T,F] representation of high and low phase of ts
HLTF <- data.frame(HL = c("H","L"), TF = c(TRUE,FALSE))

#Here the core thing happens
#If h[t+T2] < h[t] & h[t-T2] < h[t] then high tide, else low tide
h$TF <- (h$hfilt>approx(x=h$time,y=h$hfilt,xout= pmin(h$time+T2,h$time[length(h$time)]))$y)+hoffset&(h$hfilt>approx(x=h$time,y=h$hfilt,xout=pmax(h$time[1],h$time-T2))$y+hoffset)
h$HL <- HLTF$HL[match(h$TF,HLTF$TF)]

#Give every high and low tide phase a number
h$N <- 0
h$N[2:(length(h$time))] <- 1*(h$HL[1:(length(h$time)-1)] != h$HL[2:(length(h$time))])
h$N[1] <- 1
h$N <- cumsum(h$N)


#Now, find all maxima within each high and low phase
#
#Remark: a very short and clean way of coding would be like this
#
#minmax <- by(h,h$N,function(x,...){
#			switch(x$HL[1], 
#				L = x[which.min(x$h),], 
#				H = x[which.max(x$h),])},
#			simplify=T)
#HL <- do.call(rbind,minmax)
#
#However, this is about twice as slow (due to do.call() and by()) as the following, dirtier code

max <- tapply(h$h,h$N,max)
min <- tapply(h$h,h$N,min)
h$max <- max[h$N]
h$min <- min[h$N]
h$HLval <- 0		
h$HLval <- (h$max==h$h & h$HL=="H")*h$h + (h$min==h$h & h$HL=="L")*h$h   #Pick minimum in low water phase and maximum in high water phase
HL <- h[h$HLval != 0,] 	#Select only High and Low waters
HL <- HL[match(unique(HL$N),HL$N),]	#cleanup: take first occurance in cases where maximum or minimum is reached multiple times during a single water phase.


h$h <- h$ho
return(list(HL = HL[c("time","h","HL","h0")], 		# Data frame with extrema
              h = h[c("time","h","h0","HL","N")]) 	# Original water level data frame with additional attributes
              )
}


#' @name IT
#' @aliases IT
#' @title Inundation time
#' @description Calculate inundation times, i.e. time intervals for which water level h > h0. Care must be taken when there are gaps (long time periods for which there is no data )in the time series. Either the erroneous values have to removed manually, or a wrapper making use of the function gapsts can be used.
# \usage{IT(h, h0, h0marg = 0.3, dtMax, unit = "mins")}
#' @param h Water level time series. data frame with time and h column
#' @param h0 Reference level, either single valued or vector with same length as h
#' @param hoffset Offset level to cope with small fluctuations due to rain, ripples. h <= h0 + hoffset is considered dry; h> h0+hoffset is considered wet
#' @param dtMax Maximum time interval in continuous water level series. Larger time intervals are considered gaps
#' @param unit Unit of dtMax. 
#' @return  a list containing:
#' \itemize{
#'    \item{IT }{Data frame with start time (t1), end time (t2) and duration (dt, unit = unit) of inundation}
#'    \item{DT }{Data frame with start time (t1), end time (t2) and duration (dt, unit = unit) of dry time}
#' }
#' @author Tom Cox <tom.cox@uantwerp.be>
#' @keywords utilities


IT <- function(h,             #Water level time series. data frame with time and h column
                h0,           #Reference level, either single valued or vector with same length as h
                hoffset = 0,  #Margin on reference level, to cope with small fluctuations in the Water level time series
                dtMax = 15,   #Maximum time interval in continuous water level series
                unit = "mins"  #Unit of dtMax and output time intervals
                )
{

dry <- h[h$h<=(h0 + hoffset),][c("time","h")]

if (dim(dry)[1] == 0) {		
#If the site never falls dry, inundation time equals the time of the time series
  IT <- data.frame(t1 = h$time[1],t2 = h$time[length(h$time)], dt = difftime(h$time[length(h$time)],h$time[1],units = unit))
  DT <- NULL
} else 
{

wet <- h[h$h>(h0 + hoffset),][c("time","h")] 		#dry time = 'inverse' of inundation time

if (dim(wet)[1] == 0) {	
#If the site is never inundated, dry time equals the time of the time series
  IT <- NULL	
  DT <- data.frame(t1 = h$time[1],t2 = h$time[length(h$time)], dt = difftime(h$time[length(h$time)],h$time[1],units = unit))	
} else
{
if (wet$time[1] > h$time[1]) {
			wet <- rbind(h[1,],wet)
			wet$time[1] <- wet$time[1] - unclass(difftime(h$time[2],h$time[1],units="secs"))
			}
if (wet$time[length(wet$time)] < h$time[length(h$time)]){
			wet <- rbind(wet,h[length(h$time),])
#			wet$time[length(wet$time)] <- wet$time[length(wet$time)] + unclass(difftime(h$time[2],h$time[1],units="secs"))
			}
if (dry$time[1] > h$time[1]) {
			dry <- rbind(h[1,],dry)
			dry$time[1] <- dry$time[1] - unclass(difftime(h$time[2],h$time[1],units="secs"))
			}
if (dry$time[length(dry$time)] < h$time[length(h$time)]) {
			dry <- rbind(dry,h[length(h$time),])
#			dry$time[length(dry$time)] <- dry$time[length(dry$time)] + unclass(difftime(h$time[2],h$time[1],units="secs"))
			}
# Calculate Inundation Time (IT) as the gaps in the DRY time series
# Calculate Dry Time (DT) as the gaps in the WET time series
# Sum of ITs and DTs should equal total length of time series

gaps <- gapsts(h$time, dtMax, unit=unit, shiftbegin=TRUE) 
IT <- gapsts(dry$time,dtMax,unit=unit, shiftbegin=TRUE)
DT <- gapsts(wet$time,dtMax,unit=unit, shiftbegin=TRUE)

if (!is.null(gaps)){
# remove gaps that are erroneously detected as dry or wet times
for (gapi in 1:length(gaps[,1])){
  IT <- IT[!(IT$t1==gaps$t1[gapi]&IT$t2==gaps$t2[gapi]),]
  DT <- DT[!(DT$t1==gaps$t1[gapi]&DT$t2==gaps$t2[gapi]),]
}
}

# ADJUSTED
# The IT or DT (one of both) will start with a gap by default, so suppress the warnings
IT <- suppressWarnings(gapsts(dry$time,dtMax,unit=unit, shiftbegin=TRUE))
DT <- suppressWarnings(gapsts(wet$time,dtMax,unit=unit, shiftbegin=TRUE))


}
}



return(list(IT = IT,     #Data frame with start time (t1), end time (t2) and duration (dt, unit = unit) of inundation
            DT = DT))    #Data frame with start time (t1), end time (t2) and duration (dt, unit = unit) of dry time
}

IF <- function(H,                  # Vector containing high water levels. 
                h0,                # Reference level for which IF has to be calculated, either single valued or array of length = length(H[1,])
                N  = length(H[,1]) # Number of cycles in time series, equals the number of high water levels when these are complete (= default value)
                )
{
IF <- length(H[H>h0])/N
return(IF*100)               #Inundation frequence [%] at (varying) reference level h0)
}
