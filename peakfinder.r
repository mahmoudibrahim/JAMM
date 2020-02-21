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
samplingSeed = 1011414 #the seed makes the random sampling that JAMM does "deterministic". You can change this to your integer seed of choice to reproduce exact same results or randomize it this way: ceiling(rnorm(1,50,300)). 
reportNoClust = "n" #report windows for which clustering failed? Reported windows will be marked by "NoClust" in the peak name (4th column). 
presched = TRUE #Set preschedule in mclapply to TRUE or FALSE
cutoff = NA #To enforce an SNR cutoff ratio for bin enrichment calling, delete NA and enter the number you want.
strict = 1 #To make bin enrichment calling more / less strict, increase / decrease this number.
meanAdjust = "n" #Adjust the initialization mean vector for each window before clustering? If you want to do so, change this to "y".
options(warn = -1, scipen = 1000) #R will not report any warnings (warn = -1), R will not use scientific notation (scipen = 1000).
#=======================> DONE! 




# ================================ 
# Required Libraries check & load
# ================================ 
if ((is.element('mclust', installed.packages()[,1]) == FALSE) || (is.element('signal', installed.packages()[,1]) == FALSE) || (is.element('parallel', installed.packages()[,1]) == FALSE)) {
	stop("R package 'mclust', 'signal' and 'parallel' are required. Please install them!")
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
suppressPackageStartupMessages(library("mclust"))
suppressPackageStartupMessages(library("signal"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("Rsamtools"))
suppressPackageStartupMessages(library("GenomicAlignments"))
suppressPackageStartupMessages(library("GenomicRanges"))
#=======================> DONE! 





# ================= 
# Custom Functions
# =================
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




makeVectors <- function(start,end,value)
{
  return(c(rep(value,times=end-start+1)))
}




#callCountReads <- function(cornum,type,bamfiles,index,reads,frag,chromsize,filelist,chrcount,kpduplicates,RCV)
#{
#  if (type == "single") {
#    if (cornum > 1) {
#      counts = mclapply(bamfiles,countreads,index=index,reads=reads,frag=frag,chromsize=chromsize,filelist=filelist,chrcount=chrcount,mc.cores=cornum,mc.preschedule=presched,kpduplicates=kpduplicates,RCV=RCV)
#    } else {
#      counts = lapply(bamfiles,countreads,index=index,reads=reads,frag=frag,chromsize=chromsize,filelist=filelist,chrcount=chrcount,kpduplicates=kpduplicates,RCV=RCV)
#    }
#  }  
#  if (type == "paired") {
#    if (cornum > 1) {
#      counts = mclapply(bamfiles,countreadspe,index=index,reads=reads,frag=frag,chromsize=chromsize,filelist=filelist,chrcount=chrcount,mc.cores=cornum,mc.preschedule=presched,RCV=RCV)
#    } else {
#      counts = lapply(bamfiles,countreadspe,index=index,reads=reads,frag=frag,chromsize=chromsize,filelist=filelist,chrcount=chrcount,RCV=RCV)
#    }
#  }
#  return (counts)
#}




#Get per-row Geometric mean (takes list, returns vectors, not lists!)
geomeanL <- function(mat){
	n = length(mat)
	if (n > 1) {
		mult = (mat[[1]])*(mat[[2]])
		if (n > 2) {
			for (i in 3:n) {
				mult = mult*(mat[[i]])
			}
		}
		mult = mult^(1/n)
		mat = mult
		return(mat)
	} else {
		return(mat[[1]])
	}
}



#Get per-row Geometric mean (takes matrix, returns vectors, not lists!)
geomean <- function(mat){
	n = NCOL(mat)
	if (n > 1) {
		mult = (mat[,1])*(mat[,2])
		if (n > 2) {
			for (i in 3:n) {
				mult = mult*(mat[,i])
			}
		}
	mult = mult^(1/n)
	mat = mult
	}
return(mat)
}



##Read in bed(st2) file
#parsein = function(bedfile) {
#  l = read.table(bedfile, header = FALSE)[[1]]
#	l = l + 1
#	return(l)
#}


##Read in bedpe(st2) file
#parseinpe = function(bedfile) {
#	l = read.table(bedfile, header = FALSE)
#	l = cbind((l[[1]] + 1), l[[2]])
#	return(l)
#}


##Produces normalized extended read counts (takes output of parsein(), return a vector of floats)
#countreads = function(bamfile, index, reads, frag, chromsize, filelist, chrcount, kpduplicates, RCV=ReadChromVector) {
#  #o ist die Nummer des Chromosoms
#  ###########COMMENT FROM MAHMOUD################
#  ## In my implementation, o was not the chromosome, it was which replicate. Just pointing it out in case you didn't set it up here on purpose
#  ###########ANSWER##############################
#  ## well actually that was just my first guess, later I figured out what it was and that it was pretty useful, I used in all scripts to know which replicate it was
#  ###############################################
#	o = which(filelist == bamfile)
#	#print(paste("o",o,class(o)))
#  
#  cat("\n", "peakfinder: Counting",bamfile)
#  if (o > length(frag))
#  {
#    FragLength <- frag[length(frag)]
#  } else {
#    FragLength <- frag[o]
#  }
#	
#  als <- NA
#  #o = which(filelist == bamfile)
#  bamfile <- BamFile(bamfile)
#  param <- ScanBamParam()
#  #countlist <<- list()
#  countlist <- list()
#  #creates a list for all chromosomes to call peaks at
#  countlist <- makeCountList(seqnames(seqinfo(bamfile)),chromName,chromSize,RCV)
#  chromnames <- names(countlist)
#  starting=TRUE
#  curnum <- 1
#  for (element in chromnames)
#  {
#    if(curnum == 1) {
#      cat("\n","Counting:",length(chromnames),"elements","\n")
#      cat(paste0(element,", ",curnum,"; "))
#    } else {
#      cat(paste0(element,", ",curnum,"; "))
#    }
#    #print(paste("Counting",element))
#    #chromlength <- length(countlist[[element]])
#    chromlength <- as.integer(countlist[[element]])
#    
#   #####is it important to skip chromosome if length doesn't fit?
#	###########COMMENT FROM MAHMOUD################
#	## I think so, yes. Would be nice to get an error here and skip
#	###############################################
#  ## the alsover element a bit down from here would make sure that there are no reads after the chromosome end. So should I skip to the next chromosome or quit?
#  ## at the moment it's quitting but of course we could just use a next(). At least as long as it's still a for loop, it would jump to the next chromosome. 
#  ## for apply I don't know, but I guess it works there as well.
#  ###############################################
#    
#    als <- readGAlignments(bamfile,index=indexfile,param=ScanBamParam(which=GRanges(element,IRanges(1,as.integer(chromlength)))))
#  #############COMMENT FROM MAHMOUD################
#	## what is 536870912? You call the same function twice but with that difference#
#	#############ANSWER##############################
#  ## I check whether there are reads after the end of the chromosome, so I look from chromlength to 536870912. I wanted to use Max.int or something like this, because
#  ## I don't know the maximum coordinate. But for some reason this number here is the biggest that readGAlignments allows. If I go even one number higher, I get an error.
#  ## So that's because of the package used.
#  #################################################
#    alsover <- readGAlignments(bamfile,index=indexfile,param=ScanBamParam(which=GRanges(element,IRanges(as.integer(chromlength),536870912))))
#    if(!!length(alsover)){
#      #message(paste0(chromName, ", Warning: Chromosome ",element,"has no reads, Skipped!"))
#      message(paste0(element, ", Warning: Chromosome ",element,"has no reads, Skipped!"))
#      #if (chrcount == 1) {
#      if (starting) {
#        #message(paste0(chromName, ", ERROR: The first chromosome in the analysis was skipped. I can not calculate normalization factors. You can either delete this chromosome from your chromosome size file or fix the previous warning!"))
#        message(paste0(element, ", ERROR: The first chromosome in the analysis was skipped. I can not calculate normalization factors. You can either delete this chromosome from your chromosome size file or fix the previous warning!"))
#        quit(status=1)
#      } else {
#        quit()
#      }
#    }
#    rm(alsover)
#    gc()
#    if (!length(als))
#    {
#      countlist[[element]] <- NULL
#      next()
#    }
#    als <- GRanges(als)
#    #extends/truncates reads to user specified length and deletes the ones that don't match chromosome size anymore
#    als <- resize(als, FragLength)
#    #print(paste("FragLength is",FragLength))
#    starts<-start(ranges(als));ends<-end(ranges(als));seqnames<-as.character(seqnames(als));strands<-as.character(strand(als))
#    als <- mapply(removeOutOfRangeReads,starts,ends,seqnames,strands,chromlength)
#    unl <- unlist(als)
#    starts <- as.integer(unl[c(TRUE,FALSE,FALSE,FALSE)]);ends <- as.integer(unl[c(FALSE,TRUE,FALSE,FALSE)])
#    seqnames <- as.character(unl[c(FALSE,FALSE,TRUE,FALSE)]);strands <- as.character(unl[c(FALSE,FALSE,FALSE,TRUE)])
#    als <- GRanges(seqnames=Rle(seqnames),strand=Rle(strands),ranges=IRanges(start=starts,end=ends))
#    if (!kpduplicates) {
#      als <- unique(als)
#    }
#    #calculate the number of reads on each position
#    rm(starts,ends,seqnames,strands)
#    gc()
#    cov <- coverage(als)
#    curcov <- cov[element]
#    curcov <- as(curcov, "GRanges")
#    values <- score(curcov)
#    st <- start(ranges(curcov))
#    en <- end(ranges(curcov))
#    rm(curcov,cov)
#    gc()
#    ###COMMENT FROM MAHMOUD########
#    ## I couldn't follow this part right before here
#    ## Maybe you can explain to me later
#    #####################################
#    ## I have to extract start, end, ... because to remove any reads that would be out of range after resizing, I have to change
#    ## their values outside the GRanges element. Inside the element the package won't let me make any changes, with some weird explanation
#    ## like "the user isn't supposed to change this because that would destroy the sense of some attributes, ..." Yeah so I have to make a new element
#    ## The coverage function counts how many reads are at which position, for every chromosome so then I extract the current chromosome.
#    ## As GRanges it's a better object to work with, even though again it won't let me change anything, just extract. 
#    ## Values is a big vector that is compressed, because consecutive positions in the chromosome with the same
#    ## counts are summarized. And after this comment here I decompress this to a chromosome size long vector.
#  	#####################################
#    #decompress counts to a chromosome size long vector
#    reslist <- mapply(makeVectors,st,en,values,SIMPLIFY=FALSE)
#    curvector <- unlist(reslist)
#    #normalize the read counts
#    rm(reslist,st,en)
#    gc()
#    
#    
#    mCount = mean(curvector)
#    
#    #####continue here!
#    #if (chrcount == 1) {
#    if (starting) {
#      curvector = curvector/mCount
#      #####adapt for different chromosomes, as I do all at once
#      write(paste(mCount), file = paste0(out, "/norma.", o, ".info"))
#      #write(paste(mCount), file = paste0(out, "/norma.", element, ".info"))
#    } else {
#      meanCounts = mean(as.numeric(read.table(paste0(out, "/norma.", o, ".info"))[[1]]))
#      if ((mCount >  (5*meanCounts)) || (mCount <  (0.2*meanCounts))) {
#        mCount = meanCounts
#      } else {
#        write(paste(mCount), file = paste0(out, "/norma.", o, ".info"), append = TRUE)
#      }
#      curvector = curvector/mCount
#    }
#    
#    ###COMMENT FROM MAHMOUD########
#    ## like in "xcorr", we shouldn't construct this "countlist" variable. 
#    ## If I got it right, this has all the counts for the whole genome. If so, it will for sure be huge, regardless of whether you can delete the  previous "reslist", "st"...etc.
#    ## The alternative to be able to directly call 
#    #####################################
#    ## Ok so I will try to do this. Although it might be that this only solves a small part of the problem. In one run I asked about how big this 
#    ## storing object was, and it was about 1GB. "als" will be huge, probably even bigger than countlist. If I can't remove als after counting, that will
#    ## take a lot of memory as well. 
#  	#####################################
#    countlist[[element]] <- curvector
#    rm(als,curvector)
#    gc()
#    gc()
#    #als <- NA
#    starting <- FALSE
#    #curvector <- NA
#    curnum=curnum+1
#  }
#  return (countlist)
#}



#Produces normalized extended read counts (takes output of parsein(), return a vector of floats)
#countreadspe = function(bedfile, reads, chromsize, filelist, chrcount) {
#	
#	o = which(filelist == bedfile)
#	if (reads[[o]][length(reads[[o]]),2] > chromsize) {
#		message(paste0(chromName, ", Warning: Read alignments do not match chromosome length, Skipped!"))
#		
#		if (chrcount == 1) {
#			message(paste0(chromName, ", ERROR: The first chromosome in the analysis was skipped. I can not calculate normalization factors. You can either delete this chromosome from your chromosome size file or fix the previous warning!"))
#			quit()
#			system("exit 1")
#		} else {
#			quit()
#		}
#	}
#	
#	counts = vector(mode = "numeric", length = chromsize)
#	for (j in 1:length(reads[[o]][,1])) {
#		counts[(reads[[o]][j,1]):(reads[[o]][j,2])] = counts[(reads[[o]][j,1]):(reads[[o]][j,2])] + 1
#	}
#	
#	mCount = mean(counts)
#	
#	if (chrcount == 1) {
#		counts = counts/mCount
#		write(paste(mCount), file = paste0(out, "/norma.", o, ".info"))
#	} else {
#		meanCounts = mean(as.numeric(read.table(paste0(out, "/norma.", o, ".info"))[[1]]))
#		if ((mCount >  (5*meanCounts)) || (mCount <  (0.2*meanCounts))) {
#			mCount = meanCounts
#		} else {
#			write(paste(mCount), file = paste0(out, "/norma.", o, ".info"), append = TRUE)
#		}
#		counts = counts/mCount
#	}
#
#	return(counts)
#}



#find enriched bins
pickbins = function(winStart, als, binSize, chromlength, numdup, C, cutoff, strict, mCs, dCs, bkgdS) {

	if ((winStart + binSize) <= chromlength) {
		winEnd = winStart + binSize
	} else {
		winEnd = chromlength
	}
	binSizeTemp = winEnd - winStart
	tempend = winEnd - 1

	#extract subset of the background
	if (bkgdS != "None") {
		Cs = als[[numdup+1]][winStart:tempend]
		mCs = mean(Cs)
		dCs = sd(Cs)
	}
  #print(paste("winStart",winStart,"counts",Cs,"binSize",binSize,"chromSize",chromSize,"numdup",numdup,"C",C,"cutoff",cutoff,"strict",strict,"mCs",mCs,"dCs",dCs,"bkgd",bkgd))
  #if (Cs>0){
  #  print("count is >0")
  #}
	go = rep(0, numdup)
	for (g in 1:numdup) {
		mS = (mean(als[[g]][winStart:tempend])) 
		ratio = mS/dCs
    #print(paste("mS",mS,"mCs",mCs,"strict",strict,"ratio",ratio,"cutoff",cutoff))
		#print(paste("mS",mS,"mCs",mCs,"strict",strict,"ratio",ratio,"cutoff",cutoff,"winStart",winStart,"counts",Cs,"binSize",binSize,"chromlength",chromlength,"numdup",numdup,"C",C,"dCs",dCs,"bkgd",bkgd))
		#if((is.na(dCs))||(is.null(dCs))){
		#  print("dCs is null or na")
		#  print(dCs)
		#}
		#if((is.na(g))||(is.null(g))){
		#  print("g is null or na")
		#  print(g)
		#}
		#if((is.na(numdup))||(is.null(numdup))){
		#  print("numdup is null or na")
		#  print(numdup)
		#}
		#if((is.na(winStart))||(is.null(winStart))){
		#  print("winStart is null or na")
    #  print(winStart)
		#}
		#if((is.na(tempend))||(is.null(tempend))){
		#  print("tempend is null or na")
    #  print(paste("numdup",numdup))
    #  print(tempend)
		#}
    #if ((is.na(mS))||(is.null(mS))) {
    #  print("mS is null or na")
    #  print(mS)
    #  print(als[[g]][winStart:tempend])
    #  print(paste("g",g))
    #  print(paste("winStart",winStart))
    #  print(paste("tempend",tempend))
    #}
		#if ((is.na(mCs))||(is.null(mCs))) {
    #  print("mCs is null or na")
    #  print(paste("background is",bkgd))
    #  print(paste("winstart",winStart))
    #  print(paste("tempend",tempend))
    #  print(mCs)
    #  print(als[[numdup+1]][winStart:tempend])
		#}
		#if ((is.na(strict))||(is.null(strict))) {print("strict is null or na");print(strict)}
		#if ((is.na(ratio))||(is.null(ratio))) {print("ratio is null or na");print(ratio)}
		#if ((is.na(cutoff))||(is.null(cutoff))) {print("cutoff is null or na");print(cutoff)}
		if ((mS > (mCs * strict)) && (ratio > cutoff)) {
			go[g] = 1
		}
		#if (mS > (mCs * strict)){
      #print(paste("mS",mS,"mCs",mCs,"strict",strict,"ratio",ratio,"cutoff",cutoff,"winStart",winStart,"counts",Cs,"binSize",binSize,"chromlength",chromlength,"numdup",numdup,"C",C,"dCs",dCs,"bkgd",bkgd))
      #"mS 0 mCs 0 strict 1 ratio NaN cutoff 0.0376209824577413 winStart 120651 counts 0 binSize 50 chromSize 91744698 numdup 2 C  dCs 0 bkgd /tmp/tmp.nFwYfWMmlK/bkgd.26851/ext.bkgd.chrY.ctrl.bed"
		  #"mS 0 mCs 0 strict 1 ratio NaN cutoff 0.00727247426186456 winStart 108001 counts 0 binSize 250 chromlength 91744698 numdup 2 C  dCs 0 bkgd ../Jdirectorycontrol"
		
      #if(ratio > cutoff){
      #  go[g] = 1
      #}
		#}
	}
	veep = sum(go)
	return(veep)
}




#find enriched wins
pickwins = function(winStart, coffeeshopSud, als, numdup, startlist, winSizeS) {
	
	plz = which(startlist == winStart)
	winEnd = coffeeshopSud[plz]
	rWinSize = winEnd - winStart + 1


	if(rWinSize >= winSizeS) {
		mS = rep(0, numdup)
		for (g in 1:numdup) {
			mS[g] = (mean(als[[g]][winStart:winEnd])) 
		}
		veep = mean(mS)
	} else {
		veep = FALSE
	}
	
	return(veep)
}



#score windows for fast analysis
scorewindow = function(winStart, coffeeshopSud, numdup, C, bkgdS, als, startlist) {
	
	plz = which(startlist == winStart)
	winEnd = coffeeshopSud[plz]
	
	#will store peak information
	writethis = list()

	rWinSizeTemp = winEnd - winStart + 1
	
	#extract subset of the IP
	Rs = matrix(nrow = rWinSizeTemp, ncol = numdup)
	Rsr = Rs
	for (j in 1:numdup) {
		Rsr[,j] = als[[j]][winStart:winEnd]
		Rs[,j] = filtfilt(rep(1,80)/80,1,Rsr[,j])
	}
	#extract subset of the background
	if (bkgdS != "None") {
		Cs = als[[numdup+1]]
		Cmin = min(Cs[Cs > 0])
		Cs = Cs[winStart:winEnd]
		Cs = filtfilt(rep(1,80)/80,1,Cs) + Cmin #gets rid of Inf in the fold change
	} else {
		set.seed(samplingSeed)
		Cs = sample(C, rWinSizeTemp, replace = TRUE)
		Cs = filtfilt(rep(1,80)/80,1,Cs)
	}
			
	#start scoring
	signal = (geomean(Rs))
	cairo = (mean(signal)) / (mean(Cs))
	return(cairo)
}


#Initialize MClust clustering parameters
smoothcounts = function(winStart, coffeeshopSud, numdup, als, startlist) { #helper function1

	plz = which(startlist == winStart)
	winEnd = coffeeshopSud[plz]
	
	#extract subset of the IP
	Rs = matrix(0, nrow = (winEnd - winStart + 1), ncol = numdup)
	for (j in 1:numdup) {
		Rs[,j] = als[[j]][winStart:winEnd]
	}
	#smooth extended read counts
	for (j in 1:numdup) {
		Rs[,j] = filtfilt(rep(1,80)/80,1,Rs[,j])
	}	
	return(Rs)
}
cluster = function(model, sig, init, clustnummerS, noise) { #helper function2
	set.seed(samplingSeed)
	noisy = sample(noise, length(sig[,1]), replace = TRUE)
  #print(paste("noisy",head(noisy),length(noisy),class(noisy)))
  print(paste("model",head(model),length(model),class(model)))
	#print(paste("sig",head(sig),length(sig),class(sig)))
  print(paste("init",head(init),length(init),class(init)))
  print(paste("sig+noisy",head(sig+noisy),length(sig+noisy),class(sig+noisy)))
	clust = me(model, sig+noisy, init)
  print(paste("clust",head(clust[[6]])))
	bicc =  bic(model, clust$loglik, length(sig[,1]), length(sig[1,]), as.numeric(clustnummerS))
  print(paste("bicc",head(bicc)))
	out = list(bicc = bicc, param = clust$parameters)
  print(paste("out",head(out)))
	return(out)
}
initparam = function(coffeeshopNord, coffeeshopSud, numdup, als, cornumS, clustnummerS, modelnames, noise) { #main function
	n = length(coffeeshopNord)
  print(paste("n",n,"als",head(als[[1]])))
	#smooth extended read counts
	if (cornumS > 1) {
		sig = mclapply(coffeeshopNord, smoothcounts, coffeeshopSud, numdup, als, startlist = coffeeshopNord, mc.cores = cornumS, mc.preschedule = presched)
	} else {
		sig = lapply(coffeeshopNord, smoothcounts, coffeeshopSud, numdup, als, startlist = coffeeshopNord)
	}
	sig = do.call(rbind, sig) 
  print(head(sig))
  print("after sig")
	#kmeans initialization
	set.seed(samplingSeed)
	init = kmeans(sig, clustnummerS, nstart = 20)
	init = unmap(init$cluster)
  print(head(init))
  print("after init")
  print(paste("modelnames",modelnames,"sig",head(sig),"init",head(init),"clustnummerS",clustnummerS,"noise",head(noise)))
  
	if (cornumS > 1) {
		param = mclapply(modelnames, cluster, sig, init, clustnummerS, noise, mc.cores = cornumS, mc.preschedule = presched)
	} else {
		param = lapply(modelnames, cluster, sig, init, clustnummerS, noise)
	}
  print(param)
  print("after param")
	bicc = vector(mode = "numeric", length = length(modelnames))
  print("there")
	for (i in 1:length(modelnames)) {
		bicc[i] = as.numeric(param[[i]]$bicc)
	}
  print(bicc)
	bicc = which.max(bicc)
  print(numdup)
	print(bicc)
	out = list(initparam = param[[bicc]]$param, modelname = modelnames[bicc])
	return(out) 
}




#find peaks
findpeak = function(winStart, coffeeshopSud, numdup, C, param, bkgdS, resolS, als, noise, startlist, meanAdjust, clustnummerS, chromname) {
	
	#print(winStart)
	plz = which(startlist == winStart)
	winEnd = coffeeshopSud[plz]
	
	
	#will store peak information
	writethis = list()
	ccx = 1 #default is clustering didNOT work

	
	rWinSizeTemp = winEnd - winStart + 1
	
	#extract subset of the IP
	Rs = matrix(nrow = rWinSizeTemp, ncol = numdup)
	Rsr = Rs
	if (meanAdjust == "y") {
		for (j in 1:numdup) {
			Rsr[,j] = als[[j]][winStart:winEnd]
			Rs[,j] = filtfilt(rep(1,80)/80,1,Rsr[,j])
			kabel = which.max(param$init$mean[j,])
			param$initparam$mean[j,kabel] = mean(Rs[,j])
		}
	} else {
		for (j in 1:numdup) {
			Rsr[,j] = als[[j]][winStart:winEnd]
			Rs[,j] = filtfilt(rep(1,80)/80,1,Rsr[,j])
		}
	}
	
	if (resolS != "window") {
		
		#clustering (take 1)
		take = 1
		set.seed(samplingSeed)
		noisy = sample(noise, rWinSizeTemp, replace = TRUE)
		clust = em(param$modelname, Rs+noisy, param$initparam)
		clust$classification = map(clust$z)
		if (!((any(diff(clust$classification)) != 0) && (!(any(is.na(clust$classification)))))) { #clustering didn't work, take1
			
			#repeat clustering from scratch, take 2!
			set.seed(samplingSeed)
			init = kmeans(Rs, clustnummerS, nstart = 20)
			init = unmap(init$cluster)
			set.seed(samplingSeed)
			noisy = sample(noise, rWinSizeTemp, replace = TRUE)
			clust = me(param$modelname, Rs+noisy, init)
			clust$classification = map(clust$z)
			if ((any(diff(clust$classification)) != 0) && (!(any(is.na(clust$classification))))) {
				ccx = 0 #clustering worked, take2
				take = 2
			}
		} else {ccx = 0}  #clustering worked, take1
		 
		if (ccx != 1) { #clustering worked either in take1 or take2
					
			if (numdup > 1) { #check whether all components replicates agreed on clustering assignments
				cc = vector(mode = "numeric", length = numdup)
				for (g in 1:numdup) {
					cc[g] = which.max(clust$parameters$mean[g,]) #which cluster has the largest mean (this is the peak cluster, hopefully!)
				}
				ccx = sum(diff(cc))
				cluster = cc[1]
				rm(cc)
			
				if ((ccx != 0) && (take == 1)) { #not all replicates agreed? Repeat the clustering with from scratch if not already done!
					set.seed(samplingSeed)
					init = kmeans(Rs, clustnummerS, nstart = 20)
					init = unmap(init$cluster)
					set.seed(samplingSeed)
					noisy = sample(noise, rWinSizeTemp, replace = TRUE)
					clust = me(param$modelname, Rs+noisy, init)
					clust$classification = map(clust$z)
					if ((any(diff(clust$classification)) != 0) && (!(any(is.na(clust$classification))))) { #clustering worked? check whether replicates agreed take 3
						cc = vector(mode = "numeric", length = numdup)
						for (g in 1:numdup) {
							cc[g] = which.max(clust$parameters$mean[g,]) #which cluster has the largest mean (this is the peak cluster, hopefully!)
						}
						ccx = sum(diff(cc))
						cluster = cc[1]
						rm(cc)
						take = 3
					}
				}
			} else { #no replicates!
				cluster = which.max(clust$parameters$mean) #which cluster has the largest mean (this is the peak cluster, hopefully!)
			}
		}
		
		if ((ccx != 0) && (reportNoClust=="y")) { resolS = "window" } #clustering did not work and windows should be reported
			
		
		if (ccx == 0) { #clustering worked and all replicates agree on the cluster assignments
			
			#extract subset of the background
			if (bkgdS != "None") {
				Cs = als[[numdup+1]][winStart:winEnd]
				Cs = filtfilt(rep(1,80)/80,1,Cs)
			} else {
				set.seed(samplingSeed)
				Cs = sample(C, rWinSizeTemp, replace = TRUE)
				Cs = filtfilt(rep(1,80)/80,1,Cs)
			}
			
			#find region boundaries
			loc = 1:length(clust$classification)
			gmclass = cbind(loc, clust$classification)
			locPeak = gmclass[gmclass[,2] == cluster,,drop=FALSE]
			rStart = locPeak[1] #start position of the region
			rEnd = locPeak[length(locPeak[,1]),1] #end position of the region
		
			#peak resolution check
			if (resolS == "region") {
				pSize = rEnd - rStart
				signal = (geomean(Rs[rStart:rEnd,])) 
				signal2 = (signal) - (Cs[rStart:rEnd])
				gm = mean(signal2)
				summit = which.max(geomean(Rsr[rStart:rEnd,])) - 1
				will2k = wilcox.test(signal, Cs[rStart:rEnd])
				
				#Is there signal in the region above background
				if (gm > 0) {
					writethis[[1]] = rStart + winStart - 1
					writethis[[2]] = rEnd + winStart
					writethis[[3]] = paste0(chromname, ".", rStart+winStart -1)
					#writethis[[3]] = paste0(element, ".", rStart+winStart -1)
					writethis[[4]] = "1000"
					writethis[[5]] = "."
					writethis[[6]] = gm
					writethis[[7]] = will2k$p.value
					writethis[[8]] = "-1"
					writethis[[9]] = summit
				}
			} else if (resolS == "peak") {
				#find out where separate peaks are
				d = diff(locPeak[,1]) 
				d[length(d)+1] = 0		
				locPeak = cbind(locPeak, d)
				bound1 = which(locPeak[,3] > 1, arr.in=TRUE)
				bound2 = bound1 + 1
				bound = locPeak[sort(c(bound1,bound2))]
				bound = c(rStart, bound, rEnd)
				w = 1
				warum = 0
				while (w < length(bound)) {
					pStart = bound[w] + winStart - 1 
					pEnd = bound[w+1] + winStart
					pSize = pEnd - pStart
					signal = (geomean(Rs[(bound[w]):(bound[w+1]),])) 
					signal2 = (signal) - (Cs[bound[w]:bound[w+1]])
					gm = mean(signal2)
					summit = which.max(geomean(Rsr[(bound[w]):(bound[w+1]),])) - 1
					will2k = wilcox.test(signal, Cs[(bound[w]):(bound[w+1])])
										
					weil = warum * 9
          #print(paste("weil is",weil))
					#Is there signal in the region above background
					if (gm > 0) {
						writethis[[1+weil]] = pStart
						writethis[[2+weil]] = pEnd
						writethis[[3+weil]] = paste0(chromname, ".", pStart)
						#writethis[[3+weil]] = paste0(element, ".", pStart)
						writethis[[4+weil]] = "1000"
						writethis[[5+weil]] = "."
						writethis[[6+weil]] = gm
						writethis[[7+weil]] = will2k$p.value
						writethis[[8+weil]] = "-1"
						writethis[[9+weil]] = summit
					}
					w = w + 2
					warum = warum + 1			
				}
			} #peak resolution check
		} #clustering worked and all replicates agree on clustering assignments?
	} #window resolution check
	
	if (resolS == "window") {
	
		#extract subset of the background
		if (bkgdS != "None") {
			Cs = als[[numdup+1]][winStart:winEnd]
			Cs = filtfilt(rep(1,80)/80,1,Cs)
		} else {
			set.seed(samplingSeed)
			Cs = sample(C, rWinSizeTemp, replace = TRUE)
			Cs = filtfilt(rep(1,80)/80,1,Cs)
		}
		
		#calculate scores
		pSize = rWinSizeTemp
		signal = geomean(Rs) 
		signal2 = (signal) - (Cs)
		gm = mean(signal2)
		summit = which.max(geomean(Rsr)) - 1
		will2k = wilcox.test(signal, Cs)
		
		#Is there signal in the region above background
		if (gm > 0) {
			writethis[[1]] = winStart - 1
			writethis[[2]] = winEnd
			writethis[[3]] = paste0(chromname, ".", winStart -1, ".NoClust")
			#writethis[[3]] = paste0(element, ".", winStart -1, ".NoClust")
			writethis[[4]] = "1000"
			writethis[[5]] = "."
			writethis[[6]] = gm
			writethis[[7]] = will2k$p.value
			writethis[[8]] = "-1"
			writethis[[9]] = summit
		}
	} #window reporting

return(writethis)
}


#filter return value of findpeak()
processPeaks = function(peaks) {
	peaks = matrix(unlist(peaks), ncol=9, byrow=TRUE)
	peaks = peaks[peaks[,1] != FALSE,,drop=FALSE]
	peaks = data.frame(peaks)
	return(peaks)
}

#takes one chromosome, counts, calls peaks on it
peakfindermain <- function(chromname,countlist,kpduplicates=FALSE,chromreference,ifrgd,ifrgdindex,ibkgd,ibkgdindex,frags,bins)
{
  ptm <- proc.time()
  chromlength <- as.integer(countlist[[chromname]])
  curnum = which(chromreference == chromname)
  cat(paste0(chromname,", ",curnum,"; "))
  als <- NA
  #read in read information from background. If chromosome is not in background, it won't be considered for peakcalling
  if(!is.na(ibkgd)){
    als <- mclapply(ibkgd,readdata,chromname,chromlength,kpduplicates,indexfile=ibkgdindex,mc.cores=cornumS,mc.preschedule=presched,curnum=curnum)
    if(all(is.na(als))) {
      countlist[[chromname]] <- NULL
      return(NULL)
    } else {
      alsf <- mclapply(ifrgd,readdata,chromname,chromlength,kpduplicates,indexfile=ifrgdindex,mc.cores=cornumS,mc.preschedule=presched,curnum=curnum)
      if(all(is.na(alsf))){
        countlist[[chromname]] <- NULL
        message(paste("Chromosome",chromname,"is in background, but in neither of the sample files"))
        return(NULL)
      }
      #contains read information from all sample and background files
      als <- c(alsf,als)
    }
    #without background: read all sample file information for that chromosome
  } else {
    als <- mclapply(ifrgd,readdata,chromname,chromlength,kpduplicates,indexfile=ifrgdindex,mc.cores=cornumS,mc.preschedule=presched,curnum=curnum)
    if(is.na(als)) {
      countlist[[chromname]] <- NULL
      return(NULL)
    }
  }
  print("Time til reading in data")
  print(proc.time()-ptm)
  ptm <- proc.time()
  
  #count reads for every file seperately
  print(paste("length als",length(als)))
  print(paste("length frags",length(frags)))
  for(i in 1:length(als)){
    alscur<-NA;curvector<-NA
    alscur <- als[[i]]
    if(is.na(alscur)){
      next()
    }
    if (i > length(frags))
    {
      FragLength <- frags[[length(frags)]]
    } else {
      FragLength <- frags[[i]]
    }
    
    alscur <- resizetorange(als=alscur,FragLength,chromlength)
    
    #calculate the number of reads on each position
    cov <- coverage(alscur)
    curcov <- cov[chromname]
    curcov <- as(curcov, "GRanges")
    values <- score(curcov)
    st <- start(ranges(curcov))
    en <- end(ranges(curcov))
    rm(curcov,cov)
    gc()
    ###COMMENT FROM MAHMOUD########
    ## I couldn't follow this part right before here
    ## Maybe you can explain to me later
    #####################################
    ## I have to extract start, end, ... because to remove any reads that would be out of range after resizing, I have to change
    ## their values outside the GRanges element. Inside the element the package won't let me make any changes, with some weird explanation
    ## like "the user isn't supposed to change this because that would destroy the sense of some attributes, ..." Yeah so I have to make a new element
    ## The coverage function counts how many reads are at which position, for every chromosome so then I extract the current chromosome.
    ## As GRanges it's a better object to work with, even though again it won't let me change anything, just extract. 
    ## Values is a big vector that is compressed, because consecutive positions in the chromosome with the same
    ## counts are summarized. And after this comment here I decompress this to a chromosome size long vector.
    #####################################
    #decompress counts to a chromosome size long vector
    reslist <- mapply(makeVectors,st,en,values,SIMPLIFY=FALSE)
    curvector <- unlist(reslist)
    curvector <- as.numeric(curvector)
    if(length(curvector) < chromlength){
      curvector[(length(curvector)+1):chromlength] <- 0
    }
    print(paste("als no.",i,"len curvector",length(curvector)))
    print(paste("sum curvector",sum(curvector)))
    print(paste("mean curvector",mean(curvector)))
    rm(reslist,st,en)
    gc()
    #normalize counts
    mCount = mean(curvector)
    if (curnum==1 || starting) {
      curvector = curvector/mCount
      write(paste(mCount), file = paste0(outS, "/norma.", i, ".info"))
    } else {
      meanCounts = mean(as.numeric(read.table(paste0(outS, "/norma.", i, ".info"))[[1]]))
      if ((mCount >  (5*meanCounts)) || (mCount <  (0.2*meanCounts))) {
        mCount = meanCounts
      } else {
        write(paste(mCount), file = paste0(outS, "/norma.", i, ".info"), append = TRUE)
      }
      curvector = curvector/mCount
    }
    als[[i]] <- curvector
    print(paste("als no.",i,"len curvector",length(curvector)))
    print(paste("sum curvector normalized",sum(curvector)))
    print(paste("mean curvector normalized",mean(curvector)))
    rm(curvector,alscur,curvector)
    gc()
  }
  gc()
  
  #summing counts for background
  if(length(ibkgd) > 1) {
    for (i in ((length(ifrgd)+2):length(als))) {
      als[[(length(ifrgd)+1)]] <- als[[(length(ifrgd)+1)]] + als[[i]]
    }
    als[(length(ifrgd)+2):length(als)] <- NULL
  }
  print("Time til counting reads")
  print(proc.time()-ptm)
  ptm <- proc.time()
  
  # ======================= 
  # Some preliminary stuff
  # =======================
  
  #minimum peak size (only a recommendation)
  minpeak = floor(binSize / 4)
  
  #make bins vector
  print(binSize)
  bins = seq(from = 1, to = (chromlength - 1), by = binSize)
  numdup=nrepsS
  #=======================> DONE!
  
  print("Time til preliminary stuff")
  print(proc.time()-ptm)
  ptm <- proc.time()

  cat("\n","Calling peaks of chromosome",chromname)
  if (curnum == 1 || starting) {
    write(paste(samplingSeed), file = paste0(outS, "/seed.info"))
  } else {
    samplingSeed = as.numeric(read.table(paste0(outS, "/seed.info"), header = FALSE)[[1]])
  }  

  # ============================ 
  # Estimating Background Model
  # ============================ 
  if (curnum == 1 || starting){ #first chromosome, estimate bkgd (includes SNR cutoff)
    if (is.na(cutoff)) {
      if (bkgdS != "None") {
        cutoff = vector(length = numdup)
        sdC = sd(als[[numdup+1]])
        for (x in 1:numdup) {
          cutoff[x] = (mean(als[[x]]))/(sdC)
        }
        cutoff = max(cutoff)
        C = NULL
        mCs = NULL
        write(paste(c(cutoff,NA,NA)), file = paste0(outS, "/bkgd.info"), append = TRUE)
      } else {
        cutoff = vector(length = numdup)
        mmV = var(geomeanL(als))
        mmM = mean(geomeanL(als))
        sigma = log(1+((mmV) / ((mmM)^2)))
        mu = (log(mmM)) - (0.5 * (sigma))
        set.seed(samplingSeed)
        C = rlnorm(100000, mu, sqrt(sigma))
        for (x in 1:numdup) {
          cutoff[x] = (mean(als[[x]]))/(sd(C))
        }
        cutoff = max(cutoff)
        set.seed(samplingSeed)
        snow = sample(C, binSize*5, replace = TRUE)
        mCs = mean(snow)
        dCs = sd(snow)
        write(paste(c(cutoff,sigma,mu)), file = paste0(outS, "/bkgd.info"), append = TRUE)
      }
    }
  } else { #bkgd estimated from before
    bkgdInfo = read.table(paste0(outS, "/bkgd.info"), header = FALSE)
    if (is.na(cutoff)) {
      cutoff = as.numeric(bkgdInfo[[1]][1])
    }
    
    if (bkgdS != "None") {
      C = NULL
      mCs = NULL
    } else {
      sigma = as.numeric(bkgdInfo[[1]][2])
      mu = as.numeric(bkgdInfo[[1]][3])
      set.seed(samplingSeed)
      C = rlnorm(100000, mu, sqrt(sigma))
      set.seed(samplingSeed)
      snow = sample(C, binSize*5, replace = TRUE)
      mCs = mean(snow)
      dCs = sd(snow)
    }
  }
  print("Time til estimated background model")
  print(proc.time()-ptm)
  ptm <- proc.time()
  #=======================> DONE!
  
  
  ###COMMENT FROM MAHMOUD####
  ## Starting here, stuff would be groupped in one function, that just gets called directly inside the readcounting function"..that way peaks are called "on the fly" while getting reads from a file.
  ###########################
  
  # ======================== 
  # Picking Enriched Windows
  # ========================
  if (cornumS > 1) {
    coffeeshop = mclapply(bins, pickbins, als, binSize, chromlength, numdup, C, cutoff, strict, mCs, dCs, bkgdS, mc.cores = cornumS, mc.preschedule = presched)
  } else {
    coffeeshop = lapply(bins, pickbins, als, binSize, chromlength, numdup, C, cutoff, strict, mCs, dCs, bkgdS)
  }
  coffeeshop = as.numeric(unlist(coffeeshop))
  coffeeshop[coffeeshop != numdup] = 0
  
  if (sum(coffeeshop) != 0) { #Any enriched bins?
    coffeeshop = c(0, diff(coffeeshop))
    coffeeshop = cbind(coffeeshop, bins)
    coffeeshopNord = coffeeshop[coffeeshop[,1] == numdup,,drop=FALSE]
    coffeeshopSud = coffeeshop[coffeeshop[,1] == -numdup,,drop=FALSE]
    if (nrow(coffeeshopNord) != 0) { #Any enriched bins?
      coffeeshopNord = coffeeshopNord[,2]
      coffeeshopSud = coffeeshopSud[,2] - 1
      if (length(coffeeshopSud) < length(coffeeshopNord)) {
        coffeeshopSud = c(coffeeshopSud, chromlength) 
      } else if (length(coffeeshopSud) > length(coffeeshopNord)) {
        coffeeshopNord = c(1, coffeeshopNord)
      }
      if (coffeeshopSud[length(coffeeshopSud)] > chromlength) {
        coffeeshopSud[length(coffeeshopSud)] = chromlength
      }
      
      
      if (cornumS > 1) {
        coffeeshop = mclapply(coffeeshopNord, pickwins, coffeeshopSud, als, numdup, startlist = coffeeshopNord, winSizeS, mc.cores = cornumS, mc.preschedule = presched)
      } else {
        coffeeshop = lapply(coffeeshopNord, pickwins, coffeeshopSud, als, numdup, startlist = coffeeshopNord, winSizeS)
      }
      coffeeshop = as.numeric(unlist(coffeeshop))
      coffeeshop = cbind(coffeeshopNord, coffeeshopSud, coffeeshop)
      coffeeshop = coffeeshop[coffeeshop[,3] != FALSE,,drop=FALSE]
      if (nrow(coffeeshop) > 1) {
        coffeeshop = coffeeshop[order(coffeeshop[,3], decreasing = TRUE),]
      }
      rm(bins)
      print("Time til picking enriched bins")
      print(proc.time()-ptm)
      ptm <- proc.time()
      #=======================> DONE!
    
    
    
    
      # =================================== 
      # Initializing Clustering Parameters
      # ===================================
      if (length(coffeeshop[,1]) > 0) { #any enriched windows detected?
        if (initializeS == "deterministic") {
          yummy = ceiling(length(coffeeshop[,1]) / 1000)
          if (yummy == 0) {
            yummy = 1
          }
        }
        if (initializeS == "stochastic") {
          yummy = ceiling(length(coffeeshop[,1]) / 4)
          if (yummy > 20) {
            set.seed(samplingSeed)
            yummy = sample(1:yummy, 20)
          } else if (yummy > 0) {
            yummy = 1:yummy
          } else {
            yummy = 1
          }
        }
        coffeeshopNord = coffeeshop[yummy,1]
        coffeeshopSud = coffeeshop[yummy,2]
        set.seed(samplingSeed)
        noise = rnorm(100000, mean=0, sd=0.1)
        param = initparam(coffeeshopNord, coffeeshopSud, numdup, als, cornumS, clustnummerS, modelnames, noise)
        print("Time til clustering parameters")
        print(proc.time()-ptm)
        ptm <- proc.time()
        #=======================> DONE!
        
        
        
        
        
        # ========================== 
        # Enriched Window Filtering
        # ==========================
        if (windoweS != 1) { #do it only if window fold enrichment filtering is required
          if (cornumS > 1) {
            scores = mclapply(coffeeshop[,1], scorewindow, coffeeshop[,2], numdup, C, bkgdS, als, startlist = coffeeshop[,1], mc.cores = cornumS, mc.preschedule = presched)
          } else {
            scores = lapply(coffeeshop[,1], scorewindow, coffeeshop[,2], numdup, C, bkgdS, als, startlist = coffeeshop[,1])
          }
          scores = unlist(scores)
          
          if (windoweS == "auto") {
            lscores = log(scores)
            if (length(scores) > 0) {
              if (curnum == 1 || starting) {
                cutthisTEMP = ((mean(lscores)) + (sd(lscores)*1))
                write(paste(cutthisTEMP), file = paste0(outS, "/bkgd.info"), append = TRUE)
              } else {
                cutthisTEMP = as.numeric(bkgdInfo[[1]][4])
              }
              finalwins = which(lscores > cutthisTEMP)
              cutthisW = min(scores[finalwins])
              coffeeshop = cbind(coffeeshop[finalwins,1], coffeeshop[finalwins,2], coffeeshop[finalwins,3])
            } else {
              if (curnum == 1 || starting) {
                cutthisTEMP = 0
                cutthisW = "Not Applicable, All Windows Analyzed!"
                write(paste(cutthisTEMP), file = paste0(outS, "/bkgd.info"), append = TRUE)
              }
            }
          } else {
            cutthisW = windoweS
            if (length(scores) > 0) {
              finalwins = which(scores >= windoweS)
              coffeeshop = cbind(coffeeshop[finalwins,1], coffeeshop[finalwins,2], coffeeshop[finalwins,3])
            }
          }
        } else { cutthisW = 1 }
        print("Time til enriched window filtering")
        print(proc.time()-ptm)
        ptm <- proc.time()
        #=======================> DONE!
        
        
        
        
        
        # ============== 
        # Finding Peaks
        # ==============
        if (nrow(coffeeshop) > 0) { #any enriched windows left after filtering?
          coffeeshop = cbind(coffeeshop[,1], coffeeshop[,2])
          if (cornumS > 1) {
            peaks = mclapply(coffeeshop[,1], findpeak, coffeeshop[,2], numdup, C, param, bkgdS, resolS, als, noise, startlist = coffeeshop[,1], meanAdjust, clustnummerS, chromname=chromname,mc.cores = cornumS, mc.preschedule = presched)
          } else {
            peaks = lapply(coffeeshop[,1], findpeak, coffeeshop[,2], numdup, C, param, bkgdS, resolS, als, noise, startlist = coffeeshop[,1], meanAdjust, clustnummerS, chromname=chromname)
          }
          if (!(is.null(peaks))) { #any peaks discovered?
            writethis = processPeaks(peaks)
            print("Time til finding peaks")
            print(proc.time()-ptm)
            ptm <- proc.time()
            #=======================> DONE!
            
            
            
            
            
            # =========================
            # Writing Peak Information
            # =========================
          } else { nothing = TRUE } #no peaks
        } else { nothing = TRUE } #no enriched windows left after filtering
      } else { nothing = TRUE; cutthisW = windoweS } #no enriched widnows discovered
    } else { nothing = TRUE; cutthisW = windoweS } #no enriched bins discovered
  } else { nothing = TRUE; cutthisW = windoweS } #no enriched bins discovered
  
  if (isTRUE(nothing)) {
    file.create(paste0(outS, "/", chromname, ".peaks.bed"))
    print(paste0("file is: ",outS, "/", chromname, ".peaks.bed"))
    write(paste(chromname, minpeak, sep = "  "), file = paste0(outS, "/min.peaksize"), append=TRUE)
    
    if (curnum == 1 || starting) {
      message(paste0("No peaks found! - Window Fold Enrichment: ", cutthisW, " - Seed: ", samplingSeed))
    } else {
      message("No peaks found!")
    }
  } else {
    write(paste(chromname, writethis$X1, writethis$X2, writethis$X3, writethis$X4, writethis$X5, writethis$X6, writethis$X7, writethis$X8, writethis$X9, minpeak, sep = "	"), file = paste0(outS, "/", chromname, ".peaks.bed"), ncolumns = 1)
    write(paste(chromname, minpeak, sep = "	"), file = paste0(outS, "/min.peaksize"), append=TRUE)
    
    
    if (curnum == 1 || starting) {
      message(paste0("Done! - Window Fold Enrichment: ", cutthisW, " - Seed: ", samplingSeed))
    } else {
      message("Done!")
    }
  }
  starting <<- FALSE
  print("Time til writing peak information")
  print(proc.time()-ptm)
  ptm <- proc.time()
  #=======================> DONE!

}
  
resizetorange <- function(als,FragLength,chromlength){
  #extends/truncates reads to fragment length and deletes the ones that don't match chromosome size anymore
  print(FragLength)
  als <- resize(als, FragLength)
  starts<-start(ranges(als));ends<-end(ranges(als));seqnames<-as.character(seqnames(als));strands<-as.character(strand(als))
  als <- mapply(removeOutOfRangeReads,starts,ends,seqnames,strands,chromlength)
  unl <- unlist(als)
  starts <- as.integer(unl[c(TRUE,FALSE,FALSE,FALSE)]);ends <- as.integer(unl[c(FALSE,TRUE,FALSE,FALSE)])
  seqnames <- as.character(unl[c(FALSE,FALSE,TRUE,FALSE)]);strands <- as.character(unl[c(FALSE,FALSE,FALSE,TRUE)])
  als <- GRanges(seqnames=Rle(seqnames),strand=Rle(strands),ranges=IRanges(start=starts,end=ends))
  rm(starts,ends,seqnames,strands)
  gc()
  return(als)
}

#only keeps reads that don't lie behind end of chromosome or before 0
removeOutOfRangeReads <- function(start,end,seqname,strand,chromlength)
{
  if (start>0 && end<=chromlength) 
  { 
    output <- list(start,end,seqname,strand)
    return (output) 
  }
}


readdata <- function(bamfile, chromname, chromlength, kpduplicates=FALSE, indexfile, curnum)
{
  bamfile <- BamFile(bamfile)
  als <- NA; alsover <- NA
  param <- ScanBamParam()
  als <- readGAlignments(bamfile,index=indexfile,param=ScanBamParam(which=GRanges(chromname,IRanges(1,as.integer(chromlength)))))
  #############COMMENT FROM MAHMOUD################
  ## what is 536870912? You call the same function twice but with that difference
  #############ANSWER##############################
  ## I check whether there are reads after the end of the chromosome, so I look from chromlength to 536870912. I wanted to use Max.int or something like this, because
  ## I don't know the maximum coordinate. But for some reason this number here is the biggest that readGAlignments allows. If I go even one number higher, I get an error.
  ## So that's because of the package used.
  #################################################
  
  #are there reads after the chromosome end? If yes length information probably doesn't fit file
  alsover <- readGAlignments(bamfile,index=indexfile,param=ScanBamParam(which=GRanges(chromname,IRanges(as.integer(chromlength),536870912))))
  if(!!length(alsover)){
    message("Warning: Chromosome ",chromname,"has reads beyond chromosome length!")
    if (curnum==1) {
      message("ERROR: The first chromosome in the analysis was skipped. I can not calculate normalization factors. You can either delete this chromosome from your chromosome size file or fix the previous warning!")
      quit(status=1)
    } else {
      return(NA)
    }
  }
  rm(alsover)
  
  if (!length(als))
  { 
    return(NA)
  }
  
  als <- GRanges(als)
  if (kpduplicates){
    als <- unique(als)
  }
  return(als)
}

#iterates over all different chromosomes existing in all files
iterateoverchromosomes <- function (ibamS, iindexS, kpduplicates=FALSE, RCV=ReadChromVector, nbkgdS,bins)
{
  ptm <- proc.time()
  #store background and sample files and their index files in different variables
  if(nbkgdS!=0)
  {
    ibkgd <- ibamS[(length(ibamS)-nbkgdS+1):length(ibamS)]
    ibkgdindex <- iindexS[(length(iindexS)-nbkgdS+1):length(iindexS)]
    ifrgd <- ibamS[1:(length(ibamS)-nbkgdS)]
    ifrgdindex <- iindexS[1:(length(iindexS)-nbkgdS)]
    print(ibkgd);print(ifrgd);
  } else {
    ifrgd <- ibamS
    ifrgdindex <- iindexS
    ibkgd <- NA
    ibkgdindex <- NA
  }
  #collect all different chromosomes from all the files
  countlist <- list()
  for (i in 1:(length(ibamS))){
    bamfile <- BamFile(ibamS[[i]])
    countlist <- c(countlist,makeCountList(seqnames(seqinfo(bamfile)),chromName,chromSize,RCV))
    chromnames <- names(countlist)
  }
  countlist<-unique(countlist)
  chromnames <- unique(chromnames)
  names(countlist) <- chromnames
  store <- list()
  cat("\n","Counting:",length(chromnames),"elements","\n")
  print("Time til calling peakfindermain")
  print(proc.time()-ptm)
  dummy <- mclapply(chromnames,peakfindermain,countlist=countlist,chromreference=chromnames,ifrgd=ifrgd,ifrgdindex=ifrgdindex,ibkgd=ibkgd,ibkgdindex=ibkgdindex,mc.cores=cornumS,mc.preschedule=presched,frags=frags,bins=bins)
}

#=======================> DONE!









# ========================== 
# Parse-in System Variables
# ==========================
#args = commandArgs(trailingOnly = TRUE) # Read Arguments from command line

#Parsing arguments and storing values
#for (each.arg in args) {
#	#chormosome size file
#	if (grepl('-sfile=',each.arg)) {
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#				size.file <- arg.split[2]
#		} else {
#			stop('No genome size file')
#		} 
#	}
#	#bed files directory
#	if (grepl('-frag=',each.arg)) {
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#				frag <- arg.split[2]
#		} else {
#			stop('No fragment length given')
#		} 
#	}	
#	#background files directory
#	if (grepl('-bkgd=',each.arg)) {
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#				bkgd <- arg.split[2]
#		} else {
#			message("No reads found in one or more replicates!")
#			quit()
#		} 
#	}	
#	#background files directory
#	if (grepl('-out=',each.arg)) {
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#				out <- arg.split[2]
#		} else {
#			stop('No output directory given')
#		} 
#	}
#	#Cluster number
#	if (grepl('-clustnummer=',each.arg)) {
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#				clustnummer <- as.numeric(arg.split[2])
#		} 
#	}
#	#resolution
#	if (grepl('-resolution=',each.arg)) {
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#				resol <- arg.split[2]
#		} 
#	}
#	#processor cores
#	if (grepl('-p=',each.arg)) {
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#				cornum <- as.numeric(arg.split[2])
#		} 
#	}
#	#minimum window size
#	if (grepl('-window=',each.arg)) {
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#			winSize <- arg.split[2]
#		} 
#	}
#	#window size
#	if (grepl('-bin=',each.arg)) {
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#			binsize <- arg.split[2]
#		} 
#	}
#	#type (paired / single)
#	if (grepl('-type=',each.arg)) {
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#			type <- arg.split[2]
#		} 
#	}
#	#bed file names
#	if (grepl('-ibam=',each.arg)) {
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#				ibam <- arg.split[2]
#		} else {
#			message("No reads found in one or more replicates!")
#		} 
#	}
#	#index files
#	if (grepl('-iindex=',each.arg)) {
#	  arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#	  if (! is.na(arg.split[2]) ) {
#	    iindex <- arg.split[2]
#	  } else {
#	    stop('No index file for one or multiple bam files given')
#	  } 
#	}	
#	#window enrichment cutoff
#	if (grepl('-windowe=',each.arg)) {
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#			windowe <- arg.split[2]
#		} 
#	}
#	#initialize
#	if (grepl('-initModel=',each.arg)) {
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#			initialize <- arg.split[2]
#		} 
#	}
#	#keep duplicates or not
#	if (grepl('-uniq=',each.arg)) {
#	  arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#	  if (! is.na(arg.split[2]) ) {
#	    uniq <- as.character(arg.split[2])
#	    if(uniq=="y"){
#	      uniq=TRUE
#	    } else {
#	      uniq=FALSE
#	    }
#	  }
#	}
#	#number of replicates
#	if (grepl('-nreps=',each.arg)) {
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#				nreps <- as.numeric(arg.split[2])
#		} 
#	}
#  
#	#number of background files
#	if (grepl('-nbkgd=',each.arg)) {
#	  arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#	  if (! is.na(arg.split[2]) ) {
#	    nbkgd <- as.numeric(arg.split[2])
#	  }
#    else {
#      nbkgd <- 0
#    }
#	}
#}


if(!is.na(uniqS)){
  uniqS <- as.character(uniqS)
  if(uniqS=="y"){
    uniqS=TRUE
  } else {
    uniqS=FALSE
  }
}

if(is.na(sfileS)){
  stop('No chromosome size file')
}
if(is.na(ibamS)){
  stop('No sample or input files')
}
if(is.na(iindexS)){
  stop('No index file for one or multiple bam files given')
}
if(is.na(fragS)){
  stop('No fragment length given')
}
if(is.na(bkgdS)){
  stop('No background file directory')
}
if(is.na(nbkgdS)){
  nbkgdS <- 0
} else {
  nbkgdS <- as.numeric(nbkgdS)
}
if(is.na(outS)){
  stop('No output directory given')
}


clustnummerS<-as.numeric(clustnummerS)
##Parse in variables
chromosomes = read.table(sfileS, header=FALSE)
chromName = as.character(chromosomes$V1); #which chromosome
chromSize = as.numeric(chromosomes$V2); #chromosomes size
rm(chromosomes)
ReadChromVector <- chromName
#iindexS = strsplit(iindexS, ",", fixed = TRUE)[[1]]
#ibamS = as.list(strsplit(ibamS, ",", fixed = TRUE)[[1]])
#print(paste("ibam is",ibam))
numdup = length(ibamS) #number of replicates
nrepsS <- as.numeric(nrepsS)
cornumS <- as.numeric(cornumS)

if (numdup != nrepsS+nbkgdS) {
  print(paste("numdup, Anzahl der samplefiles:",numdup,"nreps:",nrepsS))
	message("No reads found in one or more replicates!")
	quit()
}


winSizeS = as.numeric(winSizeS)
binSize = as.numeric(binsizeS)
winSizeS = binSize * winSizeS
starting = TRUE

if (typeS == "single") {
	#frags = as.numeric(strsplit(fragS, ",", fixed = TRUE)[[1]])
  frags=as.numeric(fragS)
} else {
  message("Unfortunately, type paired doesn't work yet")
}

rm(bednames)

if (nrepsS > 1) {
	modelnames = c("VVV","VEV")
} else {
	modelnames = "V"
}

if (windoweS != "auto") {
	windoweS = as.numeric(windoweS)
}

nothing = FALSE #default is we found some peaks
options(stringsAsFactors = FALSE)
#=======================> DONE!








# =============== 
# Counting Reads
# ===============

iterateoverchromosomes(ibamS=ibamS, iindexS=iindexS, kpduplicates=uniqS,RCV=ReadChromVector,nbkgdS=nbkgdS,bins=bins)



