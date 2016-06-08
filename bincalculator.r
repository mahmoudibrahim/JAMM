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
shimazaki = function(bamfile, indexfile, rl, bins, maxIter, filelist, typeS) {
	#####type paired??
  ptm <- proc.time()
  o <- which(filelist==bamfile)
  bamfile <- BamFile(bamfile)
  param <- ScanBamParam()
  
  countlist <- makeCountList(seqnames(seqinfo(bamfile)),chromName,chromSize,ReadChromVector)
  countlist <- countlist[order(match(countlist,(sort(as.numeric(countlist),decreasing=TRUE))))]
  als <- list()
  iterate=0
  while((!length(als))&&(iterate<length(countlist))){
    iterate=iterate+1
    als <- readGAlignments(bamfile,index=indexfile,param=ScanBamParam(which=GRanges(names(countlist)[[iterate]],IRanges(1,as.integer(countlist[[iterate]])))))
  }
  message(paste("Largest chromosome in sample file is",names(countlist)[[iterate]]))
  
  if (!length(als)){
    stop("Bins couldn't be calculated: No reads on any chromosome in one or multiple sample files!")
  }
  als <- GRanges(als)
  ##########making a difference to the original JAMM
  #als <- unique(als)
  als <- resize(als,1)
  print("Time til reading, resizing data")
  print(proc.time()-ptm)
  ptm <- proc.time()
  readnum = length(als)
  maxChromLength <- countlist[[iterate]]
	readlen = rl[o]
	jack = o - 1
	bins = bins[(1+(jack*15)):(15+(jack*15))]

	costs = vector(mode = "numeric", length = length(bins))
	#Shimazaki procedure
	for (i in 1:length(bins)) {
		#construct the counting breaks vector
		genomevec = seq(1, maxChromLength, by = bins[i]);
		if (max(genomevec) < maxChromLength) {
			genomevec = append(genomevec, maxChromLength);
		}

		#create a vector of read counts
		if (typeS == "single") {
		  ameirah = sort(c(as.numeric(start(ranges(als)))))
		}
		if (typeS == "paired") {
			ameirah = sort(c((reads[[1]]), (reads[[2]])))
		}
		ameirah = hist(ameirah, breaks = genomevec, plot = FALSE)
		ameirah = ameirah$counts
		print(paste("bins no.",i,"len ameirah",length(ameirah)))
		print(paste("sum ameirah",sum(ameirah)))
		print(paste("mean ameirah",mean(ameirah)))
    #get cost function
		m = mean(ameirah)
		v = (sum((ameirah - m)^2)) / (length(ameirah))
		num = ((2*m) - v)
		den = ((bins[i]) * readnum)^2
		cost = -(log(abs(num)) - log(den))
		costs[i] = cost
    ameirah <- NA
    genomevec <- NA
	}
  print("Time til calculating costs, bins")
  print(proc.time()-ptm)
  ptm <- proc.time()
  als <- NA
  rm(ameirah,genomevec,als);gc()
  index = which.min(costs)
  finbin = bins[index]
	
  return(finbin)
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
  return(countlist)
}



#=======================> DONE! 



if (is.na(ibamS)){
  stop("No sample or input files")
}
if (is.na(iindexS)){
  stop("No index file for one or multiple bam files given")
}
if (is.na(sS)){
  stop("No chromosome size file")
}
if (is.na(readlS)){
  stop("No readlength calculated in xcorr")
}
if (is.na(storeFileS)){
  stop("No file to store result")
}
if (is.na(cornumS)){
  cornumS=1
}
if (is.na(fragsS)){
  stop("No Fragment lengths given")
}
if (is.na(typeS)){
  stop("No type given")
}
bins=NA

# ========================== 
# Parse-in System Variables
# ==========================
#args = commandArgs(trailingOnly = TRUE) # Read Arguments from command line


#Set arguments to default values
#ibam = NA # input bam file
#iindex = NA #input index file
#sFile = NA # chromosome size
#storeFile = NA # file to store result
#cornum = 1 # number of processors to use
#rl = NA # read length
#frags = NA # fragment lengths
#bins = NA

##Parsing arguments and storing values
#for (each.arg in args) {
#	if (grepl('^-s=',each.arg)) {			
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#			sFile <- arg.split[2] 
#		} else {
#			stop('No chromosome size file')
#		}
#	}
#	if (grepl('^-rl=',each.arg)) {			
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#			rl <- arg.split[2] 
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
#	if (grepl('^-ibam=',each.arg)) {
#	#if (grepl('^-ibed=',each.arg)) {
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#			ibam <- arg.split[2]
#		} else {
#			message('ERROR: The largest chromosome in your chromosome size file (-g) has no reads in one or more of your BAM files (-s). I can not calculate the bin size. You can either delete this chromosome from your chromosome size file or specify a bin size using -b parameter!')	
#			quit(status = 1)
#		} 
#	}
#	if (grepl('^-iindex=',each.arg)) {
#	  arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#	  if (! is.na(arg.split[2]) ) {
#	    iindex <- arg.split[2]
#	  } else {
#      stop('No index file for one or multiple bam files given')
#	  }
#	}
#	if (grepl('^-p=',each.arg)) {			
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#			cornum <- as.numeric(arg.split[2]) 
#		} else {
#			stop('No number of cores given')
#		}
#	}
#	if (grepl('^-f=',each.arg)) {			
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#			frags <- arg.split[2] 
#		} else {
#			stop('No Fragment lengths given')
#		}
#	}
#	if (grepl('^-type=',each.arg)) {			
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#			type <- arg.split[2] 
#		} else {
#			stop('No type given')
#		}
#	}
#	#number of replicates
#	if (grepl('-nreps=',each.arg)) {
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#				nreps <- as.numeric(arg.split[2])
#		} 
#	}
#	if (grepl('^-readl=',each.arg)) {
#	  arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#	  if (! is.na(arg.split[2]) ) {
#	    readl <- arg.split[2] 
#	  } else {
#	    stop('No readlength calculated in xcorr')
#	  }
#	}
#}

ptm <- proc.time()
#Read in variables
chromosomes = read.table(sS, header=FALSE)
chromSize = as.numeric(chromosomes$V2) #chromosome size
chromName = as.character(chromosomes$V1)
rm(chromosomes)
ReadChromVector <- chromName

#ibamS = strsplit(ibamS, ",", fixed = TRUE)[[1]]
#iindexS = strsplit(ibamS, ",", fixed = TRUE)[[1]]
print(readlS)
readlS = strsplit(readlS, ",", fixed = TRUE)[[1]]
print(readlS)
rlens = readlS[c(FALSE,TRUE)]
print(rlens)
rfile = readlS[c(TRUE,FALSE)]
print(rfile)
rl <- as.numeric(rlens[order(match(ibamS,rfile))])

if (length(ibamS) != nrepsS) {
	message('ERROR: The largest chromosome in your chromosome size file (-g) has no reads in one or more of your BED files (-s). I can not calculate the bin size. You can either delete this chromosome from your chromosome size file or specify a bin size using -b parameter. Exiting!')	
	quit()
}
#fragsS = as.numeric(strsplit(fragsS, ",", fixed = TRUE)[[1]])
fragsS <= as.numeric(fragsS)
#=======================> DONE! 




# ===================================================
# Shimazaki Procedure (Shimazaki and Shinomoto 2007)
# ===================================================

print(paste("Frags are",fragsS))
print(paste("rls are",rl))
for (i in 1:length(ibamS)) {
	if (fragsS[i] > rl[i]) {
		minbin = floor(fragsS[i] / 2)
		bins = c(bins, seq(minbin, minbin*15, by = minbin)) 
	} else {
		bins = c(bins, defaultBins)
	}
}
bins = bins[!is.na(bins)]
print("Time til starting shimazaki")
print(proc.time()-ptm)
#ptm <- proc.time()
bins = mclapply(ibamS, shimazaki, iindexS, rl, bins, maxIter, ibamS, typeS = typeS, mc.cores = cornumS)
bins = min(unlist(bins))
#=======================> DONE! 


# ================== 
# Write Information
# ==================
write(paste0(bins), file = paste0(storeFileS, "/binsize.txt"))
message(paste0("Binsize: ",bins))
#=======================> DONE!
