#!/usr/bin/env Rscript


# Starting from a DB of counts, we add some samples to it and then we call CNVs only on the new samples
# ./cnv_detection_ED.R -i /DB/Path/DB.file  -b bamlist.txt -o /output/path/ -v

start_time = Sys.time()

require(optparse) # http://cran.r-project.org/web/packages/optparse/optparse.pdf
require(ExomeDepth)
data(exons.hg19) # from ExomeDepth package

options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings

option_list = list(
  make_option(c("-b", "--bamlist"), action="store", default='', 
              type='character', help="Path to list of BAMs"),
  make_option(c("-s", "--singlebam"), action="store", default='',
              type='character', help="Path to a single BAM"),
  make_option(c("-i", "--inputControlData"), action="store", default='', 
              type='character', help="Path to Rdata"),
  make_option(c("-o", "--outdir"), action="store", default='./',
              type='character', help="Output directory [default %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print verbose output [default %default]"),
  make_option(c("-q", "--quiet"), action="store_false", dest="verbose",
              help="Do not print verbose output (this is the default)")
)
opt = parse_args(OptionParser(option_list=option_list))

print(opt,file=stderr())


# read list of BAMs
# to avoid writing a tryCatch I use file.exists, and set default to '', which is
# a file that never exists.
if (file.exists(opt$inputControlData)) {
    # read bam list directly into a vector (note use of $V1)
	load(opt$inputControlData)
} else {
	cat("You need to specify a valid Rdata file with -i.\n",file=stderr())
    cat(paste("The filename you specified was '",opt$inputControlData,"'.",sep=''),file=stderr())
    stop()
}



# read list of BAMs
# to avoid writing a tryCatch I use file.exists, and set default to '', which is
# a file that never exists.
if (file.exists(opt$bamlist)) {
    # read bam list directly into a vector (note use of $V1)
    bams = read.table(opt$bamlist,header=FALSE)$V1
} else if (file.exists(opt$singlebam)) {
    bams = c(opt$singlebam)
    if (opt$verbose) {
        cat(paste("Will read single BAM from: ","\n",sep=''),file=stdout())
        write.table(bams,row.names=FALSE,quote=FALSE,col.names=FALSE,file=stdout())
    }
} else {
    cat("**You need to specify a valid BAM list using -b or single BAM using -s.\n",file=stderr())
    cat(paste("The BAM list you specified was '",opt$bamlist,"'.",sep=''),file=stderr())
    cat(paste("The single BAM you specified was '",opt$singlebam,"'.",sep=''),file=stderr())
    stop()
}

# read output directory
# note right now if not specified, I stop execution. 
# an alternative is to create the dir.
# see http://stackoverflow.com/questions/4216753/check-existence-of-directory-and-create-if-doesnt-exist
if (file.exists(opt$outdir)) {
    setwd(opt$outdir)
} else {
    cat("You need to specify a valid output directory using -o.\n",file=stderr())
    cat(paste("The directory you specified was '",opt$outdir,"'.",sep=''),file=stderr())
    stop()
}

countsTarget = getBamCounts(bed.frame = exons.hg19, bam.files = bams)


if (opt$verbose) {
    cat(paste("Calculated counts\n",sep=''),file=stdout())
}


# counts is an S4 object.
# you need to cast it to a data frame (for bin length)
# AND to a matrix (for reference.count)
# and for some reason you can't cast S4 directly to matrix, only via df
countdfControls = as.data.frame(counts)
countdfTarget=as.data.frame(countsTarget)

countdfTarget$space <- NULL
countdfTarget$start <- NULL
countdfTarget$end <- NULL
countdfTarget$width <- NULL
countdfTarget$names <- NULL

countdf<-countdfControls


for (name in colnames(countdfTarget)){
	countdf[[name]]<-countdfTarget[[name]]
}

countmat<-as.matrix(countdf[,6:dim(countdf)[2]]) # remove cols 1-5 metadata

for (i in 1:dim(countmat)[2]) { 
    sample_name = colnames(countmat)[i]
    if (sample_name %in% colnames(countdfTarget)){
	    reference_list = select.reference.set(test.counts = countmat[,i], 
	        reference.count = countmat[,-i],
	        bin.length=(countdf$end-countdf$start)/1000,
	        n.bins.reduced = 10000)
	    reference_set = apply(
	        X = as.matrix(countdf[, reference_list$reference.choice]), 
	        MAR=1, FUN=sum)
	    all_exons = new('ExomeDepth', test=countmat[,i], 
	        reference=reference_set,
	        formula = 'cbind(test,reference) ~ 1')
	    all_exons = CallCNVs(x = all_exons, transition.probability=10^-4,
	        chromosome=countdf$space, start=countdf$start,
	        end=countdf$end, name=countdf$names)
	    write.table(all_exons@CNV.calls, file=paste(sample_name,".ED.cnv",sep=''), 
	        sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
	    if (opt$verbose) {
	        cat(paste("Wrote CNV calls for ",sample_name,"\n",sep=''),file=stdout())
	}
    }
}

duration = format(Sys.time() - start_time)

if(opt$verbose) {
    cat(paste("Completed execution in ",duration,"\n",sep=''),file=stdout())
}


