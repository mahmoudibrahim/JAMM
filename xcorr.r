########################################################################
# JAMMv1.0.7rev2 is a peak finder for joint analysis of NGS replicates.
# Copyright (C) 2014-2015  Mahmoud Ibrahim
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact: mahmoud.ibrahim@mdc-berlin.de
########################################################################





# ======================= 
# User-defined variables
# =======================
presched = TRUE #Set preschedule in mclapply to TRUE or FALSE
options(warn = -1, scipen = 1000) #R will not report any warnings (warn = -1), R will not use scientific notation (scipen = 1000).
#=======================> DONE! 



# ========================= 
# Required Libraries check 
# ========================= 
if ((is.element('parallel', installed.packages()[,1]) == FALSE)) {
	stop("R package'parallel' is required. Please install it!")
}
if ((is.element('Rsamtools', installed.packages()[,1]) == FALSE)) {
  stop("R package'Rsamtools' is required. Please install it!")
}
if ((is.element('GenomicRanges', installed.packages()[,1]) == FALSE)) {
  stop("R package'GenomicRanges' is required. Please install it!")
}
if ((is.element('GenomicAlignments', installed.packages()[,1]) == FALSE)) {
  stop("R package'GenomicAlignments' is required. Please install it!")
}
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("Rsamtools"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("GenomicAlignments"))
#=======================> DONE! 




# ============================================== 
# Parsing Arguments (source: phantom SPP script)
# ============================================== 
#args = commandArgs(trailingOnly = TRUE) # Read Arguments from command line


#Set arguments to default values
#ibam = NA #input bam file
#sFile = NA # chromosome size
#storeFile = NA # file to store result
ReadChromVector = NA


if(is.na(sS)){
  stop('No chromosome size file')
}
if(is.na(ibamS)){
  stop('No sample or input files')
}
if(is.na(iindexS)){
  stop('No index file for one or multiple bam files given')
}
if(is.na(bkgdS)){
  stop('No background directory')
}
if(is.na(nrepsS)){
  stop('No number of replicates given')
}
if(is.na(storeFileS)){
  stop('No file to store result')
}
if(is.na(cornumS)){
  stop('No number of cores given')
}

#for (each.arg in args) {
#	if (grepl('^-s=',each.arg)) {			
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#			sFile <- arg.split[2] 
#		} else {
#			stop('No chromosome size file')
#		}
#		chromosomes = read.table(sFile, header=FALSE)
#		chromName = as.character(chromosomes$V1) #which chromosome
#		chromSize = as.numeric(chromosomes$V2) #chromosome size
#		rm(chromosomes)
#	}
#
#	#input bam file
#	if (grepl('^-ibam=',each.arg)) {
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#			ibam <- arg.split[2]
#		} else {
#			message(paste0(chromName, ": No reads found in one or more replicates!"))
#			quit()
#		}
#	}
#	if (grepl('^-iindex=',each.arg)) {
#	  arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#	  if (! is.na(arg.split[2]) ) {
#	    iindex <- arg.split[2]
#	  } else {
#	    #iindex <- index(BamFile(ibam))
#	    stop('No index file for one or multiple bam files given')
#	    #quit()
#	  }
#	}
#	#background files directory
#	if (grepl('-bkgd=',each.arg)) {
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#				bkgd <- arg.split[2]
#		} else {
#			stop('No background file')
#		} 
#	}	
#	if (grepl('^-d=',each.arg)) {			
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#			storeFile <- arg.split[2] 
#		} else {
#			stop('No file to store result')
#		}
#	}
#	if (grepl('^-p=',each.arg)) {			
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#			cornum <- as.numeric(arg.split[2])
#		} else {
#			stop('No number of cores given')
#		}
#	}
#	#number of replicates
#	if (grepl('-nreps=',each.arg)) {
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#			nreps <- as.numeric(arg.split[2])
#		} 
#	}
#}

#chromosomes = read.table(sFile, header=FALSE)
chromosomes = read.table(sS, header=FALSE)
chromName = as.character(chromosomes$V1) #which chromosome
chromSize = as.numeric(chromosomes$V2) #chromosome size
rm(chromosomes)

#Read in variables
if(is.na(ReadChromVector))
{
  ReadChromVector <- chromName
}

#ibamS = strsplit(ibamS, ",", fixed = TRUE)[[1]]
#iindexS = strsplit(iindexS, ",", fixed = TRUE)[[1]]
numdup = length(ibamS) #number of replicates
numbkgd = numdup - nrepsS
if (bkgdS != "None") {
	nrepsS = nrepsS + 1
}

#####in case of error use default readlength in bincalculator

#=======================> DONE! 

matchChromNames <- function(actlevels,userchrom,userchromlength)
{
  userchrom <- tolower(userchrom)
  userchrom <- gsub("chr","",userchrom)
  userchrom <- gsub("seq","",userchrom)
  tempactlevels <- tolower(actlevels)
  tempactlevels <- gsub("chr","",tempactlevels)
  tempactlevels <- gsub("seq","",tempactlevels)
  temp <- match(tempactlevels,userchrom)
  if(!any(!is.na(temp)))
  {
    stop("Chromosome names in the bam file do not correspond to information in the size file")
  }
}

makeCountList <- function(actlevels, userlevels, userlevellength, binchromv)
{
  countlist <- list()
  userlevels <- tolower(userlevels)
  userlevels <- gsub("chr","",userlevels)
  userlevels <- gsub("seq","",userlevels)
  if(identical(class(binchromv),"character")) 
  {
    chromv <- tolower(binchromv)
    chromv <- gsub("chr","",chromv)
    chromv <- gsub("seq","",chromv)
    st <- match(userlevels,chromv)
    for(i in 1:length(st))
    {
      if(is.na(st[[i]])) {binchromv[[i]] <- 0} 
      else {binchromv[[i]] <- 1}
    }
  }
  tempactlevels <- tolower(actlevels)
  tempactlevels <- gsub("chr","",tempactlevels)
  tempactlevels <- gsub("seq","",tempactlevels)
  temp <- match(userlevels,tempactlevels)
  temp2 <- match(tempactlevels,userlevels)
  #false if 1)what's in the file is provided and what's provided is in the file (or more)
  if (any(is.na(temp))||any(is.na(temp2)))
  {
    actlevels <- actlevels[which(!is.na(temp2))]
    #what's in the file is provided, but there's more provided, so remove the unnecessary information
    userlevels <- userlevels[which(!is.na(temp))]
    binchromv <- binchromv[which(!is.na(temp))]
    userlevellength <- userlevellength[which(!is.na(temp))]
    temp <- Filter(Negate(function(x) is.na(unlist(x))), temp)
    if(any(is.na(temp))){
      message(paste("There are less chromosomes in the bam file than the sizes file indicates"))
    }
  }
  newbinchromv <- binchromv[order(temp)]
  newuserlevellength <- userlevellength[order(temp)]
  
  actlevels <- actlevels[which(newbinchromv==1)]
  newuserlevellength <- newuserlevellength[which(newbinchromv==1)]
  length(countlist) <- length(actlevels)
  names(countlist) <- actlevels
  for (i in 1:length(actlevels))
  {
    countlist[[actlevels[[i]]]] <- newuserlevellength[[i]]
  }
  return(countlist)
}



counthelp <- function(readstarts, breaks){
  if (!length(readstarts)){
    return(NA)
  }
  if (max(readstarts) > max(breaks)) {
    return(rep(0, length = length(breaks)))
  } else {
    return(hist(readstarts, breaks = breaks, plot = FALSE)$counts)
  }
}


normalize <- function (vector)
{
  return(vector/mean(vector))
}


xcorrelation <- function(chromname, countlist, kpduplicates=FALSE, chromreference, ifrgd, ifrgdindex,ibkgd,ibkgdindex)
{
  ptm <- proc.time()
  chromlength <- as.integer(countlist[[chromname]])
  curnum = which(chromreference == chromname)
  cat(paste0(chromname,", ",curnum,"; "))
  als <- NA
  if(!is.na(ibkgd)){
    als <- mclapply(ibkgd,readdata,chromname,chromlength,kpduplicates,indexfile=ibkgdindex,mc.cores=cornumS,mc.preschedule=presched)
    if(all(is.na(als))) {
      countlist[[chromname]] <- NULL
      return(NULL)
    } else {
      alsf <- mclapply(ifrgd,readdata,chromname,chromlength,kpduplicates,indexfile=ifrgdindex,mc.cores=cornumS,mc.preschedule=presched)
      if(all(is.na(alsf))){
        countlist[[chromname]] <- NULL
        message(paste("Chromosome",chromname,"is in background, but in neither of the sample files"))
        return(NULL)
      }
      als <- c(alsf,als)
    }
  } else {
    als <- mclapply(ifrgd,readdata,chromname,chromlength,kpduplicates,indexfile=ifrgdindex,mc.cores=cornumS,mc.preschedule=presched)
    if(all(is.na(als))) {
      countlist[[chromname]] <- NULL
      return(NULL)
    }
  }
  print("Time til reading in data")
  print(proc.time()-ptm)
  ptm <- proc.time()
  readlenchr <- NA
  if (!is.na(ibkgd)) {
    readlenchr <- mclapply(alsf,getreadleninfo,mc.cores=cornumS,mc.preschedule=presched)
  } else {
    readlenchr <- mclapply(als,getreadleninfo,mc.cores=cornumS,mc.preschedule=presched)
  }
  print("Time til getting readlength info")
  print(proc.time()-ptm)
  ptm <- proc.time()
  for(i in 1:length(als)){
    alsm<-NA;alsp<-NA;alscur<-NA;readstarts<-NA;breaks<-NA;curvector<-NA
    alscur <- als[[i]]
    if(is.na(alscur)){
      next()
    }
    alsm <- resize(alscur[[1]],1)
    alsp <- resize(alscur[[2]],1)
    readstarts <- list(start(ranges(alsp)),start(ranges(alsm)))
    breaks = seq(1, chromlength, by = 20)
    if (max(breaks) < chromlength) {
      breaks = c(breaks, chromlength);
    }
    curvector <- list()
    curvector <- mclapply(readstarts, counthelp, breaks=breaks, mc.cores=cornumS, mc.preschedule=presched)
    names(curvector) <- c("+","-")
    ################COMMENT FROM ANIKA############
    ## in this case there would be reads, but none of them would match the chromosome length, right?
    ## So would it be ok to skip this chromosome for the calculation in the sample files and ignore it for the readlengths as I currently do?
    ###############################################
    if ((all(curvector[[1]] == 0)||is.na(curvector[[1]])) && (all(curvector[[2]] == 0))||is.na(curvector[[2]])) {
      message(paste0(chromname, ", Warning: Read alignments do not match chromosome length in one or more replicates, Skipped!"))
      countlist[[chromname]]<-NULL
      return(NULL)
    }
    print(paste("als no.",i,"len curvector +,-",length(curvector[[1]]),length(curvector[[2]])))
    print(paste("sum curvector +",sum(curvector[[1]])))
    print(paste("sum curvector -",sum(curvector[[2]])))
    print(paste("mean curvector +",mean(curvector[[1]])))
    print(paste("mean curvector -",mean(curvector[[2]])))
    als[[i]] <- curvector
    rm(alsm,alsp,alscur,readstarts,curvector,breaks)
    gc()
  }
  print("Time til counting reads of all files")
  print(proc.time()-ptm)
  ptm <- proc.time()
  gc()
  
  # ================== 
  # Cross Correlation
  # ==================
  
  xcorr = list()
  if(length(ibkgd) > 1) {
    for (i in ((length(ifrgd)+2):length(als))) {
      als[[(length(ifrgd)+1)]][[1]] <- als[[(length(ifrgd)+1)]][[1]] + als[[i]][[1]]
      als[[(length(ifrgd)+1)]][[2]] <- als[[(length(ifrgd)+1)]][[2]] + als[[i]][[2]]
    }
    als[(length(ifrgd)+2):length(als)] <- NULL
  }
  
  cat("\n","Starting xcorr calculation: ",length(als)," files")
  xcorr = mclapply(als, xc, mc.cores = cornumS, mc.preschedule = presched)
  print(paste("Object als size:",format(object.size(als), units = "Mb")))
  print(paste("Object xcorr size:",format(object.size(xcorr), units = "Mb")))
  print("Time til xcorr calculation")
  print(proc.time()-ptm)
  ptm <- proc.time()
  #=======================> DONE! 
  
  
  
  # ===================== 
  # Write result to File
  # =====================
  for (i in 1:length(xcorr)) {
    bamname=ibamS[i]
    filename = strsplit(ibamS[i], "/", fixed = TRUE)[[1]]
    filename = filename[length(filename)]
    filename = substr(filename,1,nchar(filename)-4)
    if ((i == length(xcorr)) && (numbkgd != 0)) {
      filename = "ctrl"
    }
    
    correlation <- xcorr[[i]]
    if ((correlation <= 500) && (correlation >= 50)) { #only write this to file if xcorr value was plausible
      message(paste0(chromname, ", ", filename, ": Ok!"))
      print(paste0(storeFileS, "/xc.", filename, ".tab"))
      write(paste(chromname, correlation, sep = "\t"), file = paste0(storeFileS, "/xc.", filename, ".tab"), append = TRUE)
    } else {
      if (correlation) {
      } else {
        message(paste0(chromname, ", ", filename, ": Value Not Used!"))
      }
    }
  }
  print("Time til writing to file: next bincalculator")
  print(proc.time()-ptm)
  ptm <- proc.time()
  return(readlenchr)
}

getreadleninfo <- function(als){
  alsm <- als[[1]]
  alsp <- als[[2]]
  readstarts <- c(start(ranges(alsp)),start(ranges(alsm)))
  readends <- c(end(ranges(alsp)),end(ranges(alsm)))
  len <- length(readstarts)
  result <- readends-readstarts
  rm(alsm,alsp)
  gc()
  return(c(len,sum(result)))
}

readdata <- function(bamfile, chromname, chromlength, kpduplicates=FALSE, indexfile) 
{
  bamfile <- BamFile(bamfile)
  alsm <- NA; alsp <- NA
  param <- ScanBamParam()
  alsm <- readGAlignments(bamfile,index=indexfile,param=ScanBamParam(flag=scanBamFlag(isMinusStrand=TRUE),which=GRanges(chromname,IRanges(1,chromlength))))
  alsp <- readGAlignments(bamfile,index=indexfile,param=ScanBamParam(flag=scanBamFlag(isMinusStrand=FALSE),which=GRanges(chromname,IRanges(1,chromlength))))
  if ((!length(alsm)) && (!length(alsp)))
  { 
    return(NA)
  }
  if(!length(alsm)){
    alsm <- GRanges()
  } else {
    alsm <- GRanges(alsm)
  }
  if(!length(alsp)){
    alsp <- GRanges()
  } else {
    alsp <- GRanges(alsp)
  }
  if (kpduplicates){
    ##########this is making a difference to the original JAMM
    alsp <- unique(alsp)
    alsm <- unique(alsm)
  }
  return(list(alsm,alsp))
}



iterateoverchromosomes <- function (ibamS, iindexS, kpduplicates=FALSE, RCV=ReadChromVector, numbkgd)
{
  if(numbkgd!=0)
  {
    ibkgd <- ibamS[(length(ibamS)-numbkgd+1):length(ibamS)]
    ibkgdindex <- iindexS[(length(iindexS)-numbkgd+1):length(iindexS)]
    ifrgd <- ibamS[1:(length(ibamS)-numbkgd)]
    ifrgdindex <- iindexS[1:(length(iindexS)-numbkgd)]
    print(ibkgd);print(ifrgd);
  } else {
    ifrgd <- ibamS
    ifrgdindex <- iindexS
    ibkgd <- NA
    ibkgdindex <- NA
  }
  countlist <- list()
  for (i in 1:(length(ibamS))){
    bamfile <- BamFile(ibamS[i])
    countlist <- c(countlist,makeCountList(seqnames(seqinfo(bamfile)),chromName,chromSize,RCV))
    chromnames <- names(countlist)    
  }
  countlist<-unique(countlist)
  chromnames <- unique(chromnames)
  names(countlist) <- chromnames
  store <- list()
  cat("\n","xcorr: Counting:",length(chromnames),"elements","\n")
  
  readlen <- mclapply(chromnames,xcorrelation,countlist=countlist,chromreference=chromnames,ifrgd=ifrgd, ifrgdindex=ifrgdindex,ibkgd=ibkgd,ibkgdindex=ibkgdindex,mc.cores=cornumS,mc.preschedule=presched)
  length(store)<-length(ifrgd)
  for (j in 1:length(ifrgd)){
    store[[j]] <- 0
    length(store[[j]])<-2
    store[[j]][[1]]<- 0
    store[[j]][[2]]<- 0
    if (length(readlen)>0){
      for (i in 1:length(readlen)){
        if(!is.null(readlen[[i]])){
          store[[j]][[1]]<-store[[j]][[1]]+readlen[[i]][[j]][[1]]
          store[[j]][[2]]<-store[[j]][[2]]+readlen[[i]][[j]][[2]]
        }
      }
    }
  }
  rl <- mclapply(store,calculateAverageReadLength,mc.cores=cornumS,mc.preschedule=presched)
  print(rl)
  for (i in 1:length(ifrgd)){
    bamname <- ibamS[[i]]
    filename = strsplit(bamname, "/", fixed = TRUE)[[1]]
    filename = filename[length(filename)]
    filename = substr(filename,1,nchar(filename)-4)    
    print(paste0(storeFileS, "/xc.rl.", filename, ".tab"))
    write(paste(bamname, rl[[i]], sep= ","), file = paste0(storeFileS, "/xc.rl.", filename, ".tab"), append = TRUE)
  }
}



calculateAverageReadLength <- function(store)
{
  return(trunc((store[[2]]/store[[1]])+0.5))
}

xc = function(countlist) {
  if (sum(unlist(countlist)) != 0) {
    cat("\n","next")
		crossCorrelation = ccf(countlist[[2]], countlist[[1]], plot = FALSE); #xcorr
		crossCorrelation$lag = crossCorrelation$lag * 20; #correct lag for counts window
		maxCorr = which.max(crossCorrelation$acf);
		maxCorr = abs(crossCorrelation$lag[maxCorr]);
    print(maxCorr)
		return(maxCorr)
	} else {
		return(0)
	}
}


iterateoverchromosomes(ibamS=ibamS, iindexS=iindexS, kpduplicates=FALSE,RCV=ReadChromVector,numbkgd=numbkgd)



rm(xcorr);gc()
#=======================> DONE!




