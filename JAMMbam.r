#Rscript JAMMbam.r -s=/home/mibrahim/Dropbox/dev/jamm_testdata -g=/home/mibrahim/Dropbox/dev/jamm_testdata/chrSizes21.csize -o=/home/mibrahim/Dropbox/dev/jamm_testdata/bamrt265 -f=100,100 -b=200 -p=16

###MMI: Notes (22 February) - test on 2 samples, specifiying fragments and bins-- works now, deleted some commented out code

###MMI: Notes (4 June) - test on 1 sample
# -fragment length estimation works followed by binsize estimation works but then error occurs: Error in fragsS[i]/2 : non-numeric argument to binary operator Calls: source -> withVisible -> eval -> eval Execution halted
# -specifing fragment size: immediate error  Error: object 'readl' not found
# -specifying both fragment size and bin size: peak calling starts and finishes successfully- calls almost the same peaks


####--->need to at least have one file to source with all major functions









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

# this script should be run from bash

##########!



# ========================= 
# Process Input parameters
# =========================


#Defaults -- Change those if you want
mode="normal"
resol="peak"
cores="1"
window="2"
type="single"
windowe="1"
initModel="deterministic"
uniq="n"

#Defaults -- Do not change
sdir=""
gsize=""
out=""
binsize="ns"
fraglen="ns"
bkgdfile="None"
bdir="None"
ran=trunc(runif(1,0,32767))
wdir=tempdir()
#Sys.setenv(LC_ALL = C)
Sys.setenv(LANG = C)
#export LANG=C #locale defaults
#export LC_ALL=C #locale defaults



# ============================================== 
# Parsing Arguments (source: phantom SPP script)
# ============================================== 
args = commandArgs(trailingOnly = FALSE) # Read Arguments from command line
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", args[grep(file.arg.name, args)])
script.basename <- dirname(script.name)
sPath=script.basename


for (each.arg in args) {
  if (grepl('^-s=',each.arg)) {			
    arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
    if (! is.na(arg.split[2]) ) {
      sdir <- arg.split[2]
      #print(sdir)
    }
  }
  
  #genome size file
  if (grepl('^-g=',each.arg)) {
    arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
    if (! is.na(arg.split[2]) ) {
      gsize <- arg.split[2]
    }
    else {
      cat("Error: Chromosome size file not found\nExiting!\n")
      quit(status=1)
    }
  }
  if (grepl('^-o=',each.arg)) {
    arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
    if (! is.na(arg.split[2]) ) {
      out <- arg.split[2]
    }
  }
  #background files directory
  if (grepl('-c=',each.arg)) {
    arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
    if (! is.na(arg.split[2]) ) {
      bdir <- arg.split[2]
    }
  }
  if (grepl('^-m=',each.arg)) {			
    arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
    if (! is.na(arg.split[2]) ) {
      mode <- arg.split[2] 
    }
  }
  if (grepl('^-r=',each.arg)) {			
    arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
    if (! is.na(arg.split[2]) ) {
      resol <- arg.split[2]
    }
  }
  if (grepl('^-f=',each.arg)) {  		
    arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
    if (! is.na(arg.split[2]) ) {
      fraglen <- arg.split[2]
    }
  }
  if (grepl('-p=',each.arg)) {
    arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
    if (! is.na(arg.split[2]) ) {
      cores <- as.numeric(arg.split[2])
    } 
  }
  if (grepl('-w=',each.arg)) {
    arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
    if (! is.na(arg.split[2]) ) {
      window <- arg.split[2]
    } 
  }
  if (grepl('-b=',each.arg)) {
    arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
    if (! is.na(arg.split[2]) ) {
      binsize <- arg.split[2]
    } 
  }
  if (grepl('-t=',each.arg)) {
    arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
    if (! is.na(arg.split[2]) ) {
      type <- arg.split[2]
    } 
  }
  if (grepl('-e=',each.arg)) {
    arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
    if (! is.na(arg.split[2]) ) {
      windowe <- arg.split[2]
    } 
  }
  if (grepl('-i=',each.arg)) {
    arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
    if (! is.na(arg.split[2]) ) {
      initModel <- arg.split[2]
    } 
  }
  if (grepl('-d=',each.arg)) {
    arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
    if (! is.na(arg.split[2]) ) {
      uniq <- as.character(arg.split[2])
    } 
  }
  if (grepl('?=',each.arg)) {
    
  }
}


printusage <- function()
{
  cat("Welcome to JAMM v1.0.7rev2 (GNU GPLv3). Copyright (C) 2014-2015  Mahmoud Ibrahim.\n")
  cat("")
  cat("This program comes with ABSOLUTELY NO WARRANTY; for details visit http://www.gnu.org/licenses/gpl.html. This is free software, and you are welcome to redistribute it under certain conditions; visit http://www.gnu.org/licenses/gpl.html for details.\n")
  cat("")
  cat("OPTIONS:\n")
  cat("-s      directory containing Sample files (required)\n")
  cat("-g      Genome size file (required)\n")
  cat("-o      Output directory (required)\n")
  cat("-c      directory containing input or Control files\n")
  cat("-f      Fragment length(s) (default: estimated)\n")
  cat("-r      Resolution, peak or region or window (default: peak)\n")
  cat("-m      Mode, normal or narrow (default: normal)\n")
  cat("-i      clustering Initialization window selection, deterministic or stochastic (default: deterministic)\n")
  cat("-b      Bin Size (default: estimated)\n")
  cat("-w      minimum Window size (default: 2 --- Note: this means minimum_window_size = bin_size x the_value_of_-w)\n")
  cat("-e      window Enrichment cutoff, auto or any numeric value (default: 1 --- Set this to \"auto\" to estimate the window enrichment cutoff)\n")
  cat("-d      keep PCR Dupicates in single-end mode, y or n (default: n --- if -t is \"paired\", this option has no effect)\n")
  cat("-t      Type, single or paired (default: single, requires BED files. paired requires BEDPE files)\n")
  cat("-p      Number of processors used by R scripts (default: 1)\n")
  cat("")
}


if (mode=="normal") {clustno="2"}
if (mode=="narrow") {clustno="3"}

if (sdir=="" || gsize=="" || out==""){
  printusage()
  quit(status=1)
}

if (file.exists(paste0(out,"/peaks/"))){
  cat("\n\nOutput directory $out/peaks already exists. I can't override existing results!\n\n")
  quit(status=0)
}

if (fraglen == "ns" ){
  if (file.exists(paste0(out,"/xcorr/"))){
    cat ("\n\nOutput directory $out/xcorr already exists. I can't override existing results!\n\n")
    quit(status=0)
  }
}

nbkgd=NA
if (bdir != "None"){
  nbkgd = length(list.files(path=paste0(bdir,"/"),pattern="*.bam$")) #count how many background files
}
nreps = length(list.files(path=paste0(sdir,"/"),pattern="*.bam$")) #count how many sample files
#no sample files
if (nreps == 0){
  print("No sample files found!")
  quit(status=1)
}

#=======================> DONE!




# ============================= 
# Step One: Initial Processing
# =============================
cat ("\n\n============================================\nStarted JAMM Pipeline v1.0.7rev2...Hang on!\n============================================\n\n")

if (!file.exists(paste0(wdir,"/"))){
  dir.create(wdir, recursive=TRUE)
}
if (!file.exists(paste0(out,"/"))){
  dir.create(out, recursive=TRUE)
}
dir.create(paste0(wdir,"/bkgd.",ran,"/")) #directory to store background files
dir.create(paste0(wdir,"/sizes.",ran,"/")) #chromosomes and sizes
dir.create(paste0(wdir,"/samples.",ran,"/")) #store sample files

dupnum = nreps #count how many sample files


# ============================= 
# Step Two: Fragment Length
# =============================
#single-end
if ( type == "single" ){
  if (fraglen == "ns"){
    #Counting Where Reads Start and Calculating Cross Correlation
    dir.create(paste0(wdir,"/stats.",ran,"/")) #store count files
    dir.create(paste0(out,"/xcorr/")) #final xcorr results
    cat ("Calculating Fragment Length(s)...\n")
    samplelist=c()
    indexlist=c()
    sizefile=gsize

    
    #list of sample bed files and read lengths
    for (s in list.files(path=paste0(sdir,"/"),pattern="*.bam$", full.names=TRUE)){ #and for each sample file
      samplefile=basename(s)
      index=paste0(s,".bai")
      if (file.exists(index)){
        indexlist=c(indexlist,index)
      } else {
        cat("\n\nPlease provide index file ",index,"!\n\n",sep="")
        quit(status=0)
      }
       samplelist=c(samplelist,paste0(sdir,"/",samplefile))
    }
    
    #list of control bed files and read lengths
    if (bdir != "None" ){
      for (s in list.files(path=paste0(bdir,"/"),pattern="*.bam$", full.names=TRUE)){ #and for each sample file
        samplefile=basename(s)
        index=paste0(s,".bai")
        if (file.exists(index)){
          indexlist=c(indexlist,index)

        } else {
          cat("\n\nPlease provide index file ",index,"!\n\n")
          quit(status=0)
        }
        samplelist=c(samplelist,paste0(bdir,"/",samplefile))
      }
    }
    sS=gsize; ibamS=samplelist; iindexS=indexlist; storeFileS=paste0(wdir,"/stats.",ran); cornumS=cores; nrepsS=nreps; bkgdS=bdir
    source(paste0(sPath,"/xcorr.r"))

    #report xcorr results (samples)
    readl=""
    for (f in list.files(path=paste0(sdir,"/"),pattern="*.bam$", full.names=TRUE)){
      file=basename (f)
      samplefile = substr(file, 1, nchar(file)-4)
      dir.create(paste0(out,"/xcorr/",samplefile,"/")) #final xcorr results
      if (file.exists(path=paste0(wdir,"/stats.",ran,"/xc.",samplefile,".tab"))){
        file.copy(paste0(wdir,"/stats.",ran,"/xc.",samplefile,".tab"), paste0(out,"/xcorr/",samplefile,"/shifts.txt"))
      }
      readlfile=paste0(wdir,"/stats.",ran,"/xc.rl.",samplefile,".tab")
      
      if (file.exists(readlfile)){
        conn <- file(readlfile, open="r")
        lines <-readLines(conn)
        for (i in 1:length(lines)){
          readl=paste0(readl,",",lines[[i]])
        }
        close(conn)
      }
      infileS <- paste0(out,"/xcorr/",samplefile,"/shifts.txt"); outS <- paste0(out,"/xcorr/",samplefile)
      source(paste0(sPath,"/xcorrhelper.r"))
    }
    readl=substr(readl, 2, nchar(readl))
    if (bdir != "None"){
      dir.create(paste0(out,"/xcorr/ctrl/")) #final xcorr results
      if (file.exists(paste0(wdir,"/stats.",ran,"/xc.ctrl.tab"))){
        file.copy(paste0(wdir,"/stats.",ran,"/xc.ctrl.tab"), paste0(out,"/xcorr/ctrl/shifts.txt"))
      }
      infileS=paste0(out,"/xcorr/ctrl/shifts.txt"); outS=paste0(out,"/xcorr/ctrl")
      source(paste0(sPath,"/xcorrhelper.r"))
    }
  }
}


############paired-end----right now doesn't work?#######################################
if (type == "paired"){
  
  cat("Getting Average Fragment Length(s)...\n")
  dir.create(paste0(out,"/xcorr/")) #final xcorr results
  for (f in list.files(path=paste0(sdir,"/"),pattern="*.bed$", full.names=TRUE)){
    file=basename (f)
    samplefile = substr(file, 1, nchar(file)-4)
    dir.create(paste0(out,"/xcorr/",samplefile))
   
  }
  if (bdir != "None"){
    dir.create(paste0(out,"/xcorr/ctrl")) #final xcorr results
    for (f in list.files(path=paste0(wdir,"/bkgd.",ran,"/"),pattern="ctrl.bed$", full.names=TRUE)){ #and for each sample file
      file=basename (f)
      samplefile= substr(file, 1, nchar(file)-4)
        }
  }
}
#=======================> DONE!




# ================================= 
# Step Three: Calculating Bin Size
# =================================
if (binsize == "ns"){
  print ("Getting Bin Size: ")  
  samplelist=c()
  indexlist=c()
  frag=c()
  if (fraglen != "ns"){
    frag=fraglen
    k=1
  }
  #list of sample bed files and read lengths
  for (s in list.files(path=paste0(sdir,"/"),pattern="*.bam$", full.names=TRUE)){ #and for each sample file
    samplefile=basename (s)
    index=paste0(s,".bai")
    if (file.exists(index)){
      indexlist=c(indexlist,index)
    } else {
      cat("\n\nPlease provide index file $index!\n\n")
      quit(status=0)
    }
    samplename= substr(samplefile, 1, nchar(samplefile)-4)
    samplelist=c(samplelist,paste0(sdir,"/",samplefile))
    
    if (fraglen == "ns"){
      conn <- file(paste0(out,"/xcorr/",samplename,"/xcorrsummary.txt"), open="r")
      lines <-readLines(conn)
      for (i in 1:length(lines)){
        if(length(grep("Fragment Length:", lines[i]))>0){
          shift= substr(lines[[i]],17,nchar(lines[[i]]))
          frag=c(frag,shift)
        }
      }
      close(conn)
    }
  }
  ibamS=samplelist; iindexS=indexlist; sS=gsize; storeFileS=wdir; cornumS=cores; fragsS=frag; typeS=type; nrepsS=nreps; readlS=readl ########!!!!!!!!!!!!here is the readl issue!--> need to find out how and where read lengths are calculated
  source(paste0(sPath,"/bincalculator.r"))
}

if (binsize != "ns"){
  cat("You set a Bin Size:",binsize,"\n")
}
#=======================> DONE!



# =========================== 
# Step Four: Calling Peaks
# ===========================
dir.create(paste0(wdir,"/peaks.",ran,"/")) #store count files
dir.create(paste0(out,"/peaks/")) #store peak files

##########try with more than one sample file

cat ("Calling Peaks...(mode: ",mode,", resolution: ",resol,")\n",sep="")


#single-end reads
if (type == "single"){
  if (binsize == "ns"){
    binsize= read.table(paste0(wdir,"/binsize.txt"))[[1]]
    cat("binsize=", binsize)
  }
  counting=1;  		
  samplelist=c()
  indexlist=c()
  frag=c()
  k=1
  if (fraglen != "ns"){
    frag=fraglen
  }
  
  sizefile=gsize
  
  
  #list of sample bed files and fragment lengths
  for (s in list.files(path=paste0(sdir,"/"),pattern="*.bam$", full.names=TRUE)){ #and for each sample file
    samplefile=basename (s)
    index=paste0(s,".bai")
    if (file.exists(index)){
      indexlist=c(indexlist,index)
    } else {
      cat("\n\nPlease provide index file ",index,"!\n\n", sep="")
      quit(status=0)
    }
    
    samplelist=c(samplelist,paste0(sdir,"/",samplefile))
    samplename= substr(samplefile, 1, nchar(samplefile)-4)
    
    if (fraglen == "ns"){
      conn <- file(paste0(out,"/xcorr/",samplename,"/xcorrsummary.txt"), open="r")
      lines <-readLines(conn)
      for (i in 1:length(lines)){
        if(length(grep("Fragment Length:", lines[[i]]))>0){
          shift= substr(lines[[i]],17,nchar(lines[[i]]))
          frag=c(frag,shift)
        }
      }
      close(conn)
      k=k+1
    }
    
    k=1
    if (fraglen != "ns"){
      shift = as.numeric(unlist(strsplit(frag, ",", fixed = TRUE))[k])
      k=k+1
    }
    print(shift)
   }
  
  
  #control file
  if (bdir != "None"){
    if (fraglen == "ns"){
      conn <- file(paste0(out,"/xcorr/ctrl/xcorrsummary.txt"), open="r")
      lines <-readLines(conn)
      for (i in 1:length(lines)){
        if(length(grep("Fragment Length:", lines[[i]]))>0){
          bshift= substr(lines[[i]],17,nchar(lines[[i]]))
          frag=c(frag,bshift)
        }
      }
      close(conn)
    }
    if (fraglen != "ns"){
      l=dupnum+1
      bshift = strsplit(frag, ",", fixed = TRUE)[[l]]
    }
   
    for (s in list.files(path=paste0(bdir,"/"),pattern="*.bam$", full.names=TRUE)){ #and for each sample file
      samplefile=basename (s)
      index=paste0(s,".bai")
      if (file.exists(index)){
        indexlist=c(indexlist,index)
      } else {
        cat("\n\nPlease provide index file ",index,"!\n\n")
        quit(status=0)
      }
      samplelist=c(samplelist,paste0(bdir,"/",samplefile))
    }
  }
  #remove leading comma
  print(paste("frag is",frag))
  print(paste("uniq is",uniq))
  print(paste("out is",out))
  print(paste("clustno is",clustno))
  print(paste("resol is",resol))
  print(paste("window is",window))
  print(paste("cores is",cores))
  print(paste("binsize is",binsize))
  print(paste("type is",type))
  print(paste("initModel is",initModel))
  print(paste("windowe is",windowe))
  print(paste("nbkgd is",nbkgd))
  
  #call the peak calling R script
  sfileS=gsize; ibamS=samplelist; fragS=shift; bkgdS=bdir; outS=paste0(wdir,"/peaks.",ran,"/"); clustnummerS=clustno; resolS=resol; winSizeS=window; cornumS=cores; binsizeS=binsize; typeS=type; initializeS=initModel; windoweS=windowe; nrepsS=nreps; nbkgdS=nbkgd; uniqS=uniq; iindexS=indexlist
  source(paste0(sPath,"/peakfinder.r"))

  
  for (s in list.files(path=paste0(wdir,"/peaks.",ran,"/"),pattern="*.peaks.bed$", full.names=TRUE)){
    file=basename (s)
    file.copy(s, paste0(out,"/peaks/",file))
    file.remove(s)
  }

  counting=1;                     
}
                                                                
#paired-end reads [not working for now?]
#since the whole paired-end doesn't work, there are many parts in it that can't possibly work since we now read in bam files
#this would have to be adapted similar to single-end
if (type == "paired"){
  if (binsize == "ns"){
    conn <- file(paste0(wdir,"/binsize.txt"), open="r")
    binsize <-readLines(conn)
    close(conn)
  }
  counting=1;
  for (f in list.files(path=paste0(wdir,"/sizes.",ran,"/"), full.names=TRUE)){#,pattern="*")) #for each chromosome
    samplelist=c()
    sizefile=basename(f)
  
    #list of sample bed files and fragment lengths
    for (s in list.files(path=paste0(wdir,"/samples.",ran,"/"),pattern="*.bed$", full.names=TRUE)){ #and for each sample file
      samplefile=basename(s)

      samplelist=c(samplelist,paste0(wdir,"/samples.",ran,"/",samplefile))
      samplename=samplename= substr(samplefile, 1, nchar(samplefile)-4)

      x=paste0(sdir,"/",samplename)
      ##fi
    }
    #control file
    bkgdfile="None"
    if (bdir != "None"){
      bkgdfile=paste0(wdir,"/bkgd.",ran,"/bkgd.ctrl.bed")
    }
   
    sfileS=gsize; ibamS=samplelist; fragS=frag; bkgdS=bdir; outS=paste0(wdir,"/peaks.",ran,"/"); clustnummerS=clustno; resolS=resol; winSizeS=window; cornumS=cores; binsizeS=binsize; typeS=type; initializeS=initModel; windoweS=windowe; nrepsS=nreps; nbkgdS=nbkgd; uniqS=uniq; iindexS=indexlist
    source(paste0(sPath,"/peakfinder.r"))
    counting=counting+1

    if (file.exists(paste0(wdir,"/peaks.",ran,"/peaks.bed"))){
      file.copy(paste0(wdir,"/peaks.",ran,"/peaks.bed"),paste0(out,"/peaks/peaks.bed"))
      file.remove(paste0(wdir,"/peaks.",ran,"/peaks.bed"))
    }
  }
  counting=1
}


#concatenate, sort and filter
for (s in list.files(path=paste0(out,"/peaks/"),pattern="*.bed$", full.names=TRUE)){
  if(file.info(s)$size == 0){
    next()
  } else {
    conn <- file(s, open="r")
    write(readLines(conn), file=paste0(out,"/peaks/all.narrowPeak"), ncolumns=1, append=TRUE)
    close(conn)
  }
}


if (file.exists(paste0(out,"/peaks/all.narrowPeak")) && file.info(paste0(out,"/peaks/all.narrowPeak"))$size != 0){

  file.copy(paste0(wdir,"/peaks.",ran,"/min.peaksize"), paste0(out,"/peaks/min.peaksize"))
  bednamesS=paste0(out,"/peaks/all.narrowPeak")
  source (paste0(sPath,"/peakhelper.r"))
  bednamesS=paste0(out,"/peaks/all.narrowPeak"); destfileS=paste0(out,"/peaks/filtered.peaks.narrowPeak"); abovezerofileS=paste0(out,"/peaks/all.peaks.narrowPeak")
  source(paste0(sPath,"/peakfilter.r"))

  file.remove(paste0(out,"/peaks/all.narrowPeak"))
  file.remove(paste0(out,"/peaks/min.peaksize"))
  for (s in list.files(path=paste0(out,"/peaks/"),pattern="*.bed$", full.names=TRUE)){
    file.remove(s)
  }
} else {
  cat("\nNo Peaks Found!")
  dummy <- file.create(paste0(out,"/peaks/all.narrowPeak"))
}
#=======================> DONE!


unlink(wdir, recursive=TRUE)




cat ("\n\n========================================\nWe're done...Congratulations!\n========================================\n\n")
                                                                                                                     
