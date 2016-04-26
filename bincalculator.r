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
defaultBins = seq(50, 50*15, by = 50) # default binsize search space, used when fragment length is equal or less than read length 
#=======================> DONE! 




# ================================ 
# Required Libraries check & load
# ================================ 
if ((is.element('parallel', installed.packages()[,1]) == FALSE)) {
	stop("R package'parallel' is required. Please install it!")
}
if ((is.element('Rsamtools', installed.packages()[,1]) == FALSE)) {
  stop("R package'Rsamtools' is required. Please install it!")
}
if ((is.element('GenomicAlignments', installed.packages()[,1]) == FALSE)) {
  stop("R package'GenomicAlignments' is required. Please install it!")
}
if ((is.element('GenomicRanges', installed.packages()[,1]) == FALSE)) {
  stop("R package'GenomicRanges' is required. Please install it!")
}
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("Rsamtools"))
suppressPackageStartupMessages(library("GenomicAlignments"))
suppressPackageStartupMessages(library("GenomicRanges"))
#=======================> DONE! 




# ================= 
# Custom Functions
# =================
#Implements the Shimazaki procedure
#shimazaki = function(bamfile, indexfile, rl, bins, maxIter, filelist, maxChr, maxSize, type) {
shimazaki = function(bamfile, indexfile, rl, bins, maxIter, filelist, type) {
	#####type paired??
  o <- which(filelist==bamfile)
  bamfile <- BamFile(bamfile)
  param <- ScanBamParam()
  #if(maxChr[o]!=chromName)
  #{
  #  warning(paste("The biggest chromosome in the sizefile has no reads, using the biggest chromosome in the bam file:"),maxChr[o])
  #}
  #print(bamfile);print(chromSize[o]);print(chromName[o])
  #alsm <- readGAlignments(bamfile,index=indexfile,param=ScanBamParam(flag=scanBamFlag(isMinusStrand=TRUE),which=GRanges(maxChr[o],IRanges(1,as.integer(maxSize[o])))))
  #alsp <- readGAlignments(bamfile,index=indexfile,param=ScanBamParam(flag=scanBamFlag(isMinusStrand=FALSE),which=GRanges(maxChr[o],IRanges(1,as.integer(maxSize[o])))))
  #als <- readGAlignments(bamfile,index=indexfile,param=ScanBamParam(which=GRanges(maxChr[o],IRanges(1,as.integer(maxSize[o])))))
  als <- readGAlignments(bamfile,index=indexfile,param=ScanBamParam(which=GRanges(chromName,IRanges(1,as.integer(chromSize)))))
  
  #print(paste("maxSize[o]",maxSize[o],"maxChr[o]",maxChr[o]))
  
  #if ((!length(alsm)) && (!length(alsp)))
  if (!length(als))
  { 
    stop(paste("Bins couldn't be calculated: the largest chromosome",chromName,"has no reads!"))
  }  
  #alsm <- GRanges(alsm)
  #alsp <- GRanges(alsp)
  als <- GRanges(als)
  #if (rmduplicates) {
  #  alsm <- unique(alsm)
  #  alsp <- unique(alsp)
  #}
  als <- unique(als)
  #alsm <- resize(alsm,1)
  als <- resize(als,1)  
	#readnum = length(alsp)+length(alsm)
  readnum = length(als)
  #print(paste("readnum",readnum))
	#o = which(filelist == bamfile)
	readlen = rl[o]
	jack = o - 1
	bins = bins[(1+(jack*15)):(15+(jack*15))]

	costs = vector(mode = "numeric", length = length(bins))
	#Shimazaki procedure
	for (i in 1:length(bins)) {
		#construct the counting breaks vector
		genomevec = seq(1, chromSize, by = bins[i]);
		if (max(genomevec) < chromSize) {
			genomevec = append(genomevec, chromSize);
		}

		#create a vector of read counts
		if (type == "single") {
			#ameirah = sort(c((as.numeric(reads[reads[,2] == "+",,drop = FALSE][,1])), ((as.numeric(reads[reads[,2] == "-",,drop = FALSE][,1])) + readlen - 1)))
		  #ameirah = sort(c((as.numeric(start(ranges(alsp)))), ((as.numeric(start(ranges(alsm)))))))
		  ameirah = sort(c(as.numeric(start(ranges(als)))))
		}
		if (type == "paired") {
			ameirah = sort(c((reads[[1]]), (reads[[2]])))
		}
		ameirah = hist(ameirah, breaks = genomevec, plot = FALSE)
		ameirah = ameirah$counts
		#get cost function
		m = mean(ameirah)
		v = (sum((ameirah - m)^2)) / (length(ameirah))
		num = ((2*m) - v)
		den = ((bins[i]) * readnum)^2
		cost = -(log(abs(num)) - log(den))
		costs[i] = cost
    #alsm <- NA; alsp <- NA
    ameirah <- NA
    genomevec <- NA
	}
  
  #alsm <- NA; alsp <- NA
  als <- NA
  rm(ameirah,genomevec,als);gc()
  index = which.min(costs)
  finbin = bins[index]
	
  return(finbin)
}

#####maxchr: find a way to make that work fast
#processChromosomesInFile <- function(bamfile, indexfile, chromosomeinfo, mc.cores)
#{
#  print("Calculating average readlength")
#  countlist <- list()
#  bamfile <- BamFile(bamfile)
#  countlist <- makeCountList(seqnames(seqinfo(bamfile)),as.character(chromosomeinfo$V1),as.integer(chromosomeinfo$V2),as.character(chromosomeinfo$V1))
#  chromnames <- names(countlist)
#  print(chromnames);print(bamfile);print(indexfile)
#  als <- NA
#  store <- 0
#  for (element in chromnames)
#  {
#    param <- ScanBamParam()
#    chromlength <- as.integer(countlist[[element]])
#    print(chromlength)
#    als <- readGAlignments(bamfile,index=indexfile,param=ScanBamParam(which=GRanges(element,IRanges(1,chromlength))))
#    if (!length(als))
#    {
#      countlist[[element]] <- NULL
#      next()
#    }
#    als <- GRanges(als)
#    starts <- start(ranges(als)); ends <- end(ranges(als))
#    len <- length(starts)
#    result <- ends-starts
#    #countlist[[element]] <<- c(len, sum(result))
#    store <- store+c(len,sum(result))
#    als <- NA
#  }
#  rl <- calculateAverageReadLength(store)
#  maxc <- calculateMaxChrom(countlist)
#  countlist <- NA
#  rm(countlist);gc()
#  #return(round(store[[2]]/store[[1]]))
#  return(list(rl,maxc))
#}

#calculateAverageReadLength <- function(store)
#{
#  return(trunc((store[[2]]/store[[1]])+0.5))
#}

#calculateMaxChrom <- function(countlist)
#{
#  index <- which.max(countlist)
#  maxName <- names(index)
#  maxSize <- countlist[[index]]
#  return(list(maxName,maxSize))
#}

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
      warning("There are less chromosomes in the bam file than the sizes file indicates")
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
  #countlist <- mapply(defineLength,countlist,actlevels,newuserlevellength)
  return(countlist)
}

defineLength <- function(countlist,chromname,chromlength)
{
  countlist[[chromname]] <- chromlength
  return(countlist)
}

#=======================> DONE! 




# ========================== 
# Parse-in System Variables
# ==========================
args = commandArgs(trailingOnly = TRUE) # Read Arguments from command line


#Set arguments to default values
ibam = NA # input bam file
iindex = NA #input index file
sFile = NA # chromosome size
storeFile = NA # file to store result
cornum = 1 # number of processors to use
rl = NA # read length
frags = NA # fragment lengths
bins = NA

#Parsing arguments and storing values
for (each.arg in args) {
	if (grepl('^-s=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			sFile <- arg.split[2] 
		} else {
			stop('No chromosome size file')
		}
	}
	if (grepl('^-rl=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			rl <- arg.split[2] 
		}
	}
	if (grepl('^-d=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			storeFile <- arg.split[2] 
		} else {
			stop('No file to store result')
		}
	}
	if (grepl('^-ibam=',each.arg)) {
	#if (grepl('^-ibed=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			ibam <- arg.split[2]
		} else {
			message('ERROR: The largest chromosome in your chromosome size file (-g) has no reads in one or more of your BAM files (-s). I can not calculate the bin size. You can either delete this chromosome from your chromosome size file or specify a bin size using -b parameter!')	
			quit(status = 1)
		} 
	}
	if (grepl('^-iindex=',each.arg)) {
	  arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
	  if (! is.na(arg.split[2]) ) {
	    iindex <- arg.split[2]
	  } else {
	    #iindex <- index(BamFile(ibam))
      stop('No index file for one or multiple bam files given')
	    #quit()
	  }
	}
	if (grepl('^-p=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			cornum <- as.numeric(arg.split[2]) 
		} else {
			stop('No number of cores given')
		}
	}
	if (grepl('^-f=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			frags <- arg.split[2] 
		} else {
			stop('No Fragment lengths given')
		}
	}
	if (grepl('^-type=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			type <- arg.split[2] 
		} else {
			stop('No type given')
		}
	}
	#number of replicates
	if (grepl('-nreps=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				nreps <- as.numeric(arg.split[2])
		} 
	}
	if (grepl('^-readl=',each.arg)) {
	  arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
	  if (! is.na(arg.split[2]) ) {
	    readl <- arg.split[2] 
	  } else {
	    stop('No readlength calculated in xcorr')
	  }
	}
}

#Read in variables
chromosomes = read.table(sFile, header=FALSE)
chromSize = as.numeric(chromosomes$V2) #chromosome size
chromSize = max(chromSize) #get maximum chrom size
chromName = as.character(chromosomes$V1[[which.max(chromSize)]])
#rm(chromosomes)
#####
ReadChromVector <- chromName

ibam = strsplit(ibam, ",", fixed = TRUE)[[1]]
iindex = strsplit(ibam, ",", fixed = TRUE)[[1]]
readl = strsplit(readl, ",", fixed = TRUE)[[1]]
#print(readl)
rlens = readl[c(FALSE,TRUE)]
rfile = readl[c(TRUE,FALSE)]
rl <- rlens[order(match(ibam,rfile))]
#print(paste("ibam",ibam,"readl",readl,"rlens",rlens,"rfile",rfile,"rl",rl))

if (length(ibam) != nreps) {
	message('ERROR: The largest chromosome in your chromosome size file (-g) has no reads in one or more of your BED files (-s). I can not calculate the bin size. You can either delete this chromosome from your chromosome size file or specify a bin size using -b parameter. Exiting!')	
	quit()
}
#print(paste("rl",rl))
#rl = as.numeric(strsplit(rl, ",", fixed = TRUE)[[1]])
frags = as.numeric(strsplit(frags, ",", fixed = TRUE)[[1]])
#=======================> DONE! 




# ===================================================
# Shimazaki Procedure (Shimazaki and Shinomoto 2007)
# ===================================================


#results <- unlist(mclapply(ibam,processChromosomesInFile,indexfile=iindex,chromosomeinfo=chromosomes,mc.cores=cornum))
#rl <- results[c(TRUE,FALSE,FALSE)];maxChr<-results[c(FALSE,TRUE,FALSE)];maxSize<-results[c(FALSE,FALSE,TRUE)]
#print(paste("Average readlengths are",rl))
#print(paste("Biggest chromosome is",maxChr))
for (i in 1:length(ibam)) {
  #rl <- calculateAverageRL(ibam[[i]],iindex[[i]],chromosomes)
	if (frags[i] > rl[i]) {
		minbin = floor(frags[i] / 2)
		bins = c(bins, seq(minbin, minbin*15, by = minbin)) 
	} else {
		bins = c(bins, defaultBins)
	}
}
bins = bins[!is.na(bins)]
#print(paste("bins are",bins))
#bins = mclapply(ibam, shimazaki, iindex, rl, bins, maxIter, ibam, maxChr, maxSize, type = type, mc.cores = cornum)
bins = mclapply(ibam, shimazaki, iindex, rl, bins, maxIter, ibam, type = type, mc.cores = cornum)
bins = min(unlist(bins))
#=======================> DONE! 


# ================== 
# Write Information
# ==================
write(paste0(bins), file = paste0(storeFile, "/binsize.txt"))
#print(paste("storefile is",storeFile,"/binsize.txt"))
#message(bins)
#message(paste0("Binsize: ",bins,", largest chromosome ",maxChr))
message(paste0("Binsize: ",bins,", largest chromosome ",chromName))
#=======================> DONE!
