SKIP<-5 ##number of lines to skip to get the information
##get the number of rows to read in
filename<- "gencode.v16.annotation.gtf"
prefix<-"gencode.v16.annotation"
nrow=system(paste("wc -l",filename),intern=TRUE)
nrow=as.numeric(strsplit(nrow," ")[[1]][2])-SKIP

##read in table
tab <-  read.delim(filename,skip=SKIP,header=FALSE,stringsAsFactors=FALSE,colClasses=c("character","character","character","integer","integer","character","character","character"),comment.char="",nrow=nrow,quote="")

##split by gene, transcript, exon
Indexes <- split(1:nrow(tab),tab[,3])
geneIndex <- Indexes[["gene"]]
txIndex <- Indexes[["transcript"]]
##not using exon yet
##exonIndex <- Indexes[["exon"]]

###Make gene table and get a unique ID
tmptab <- tab[geneIndex,]
###get a unique ID
id <- strsplit(tmptab[,9],";")
id <- sapply(id,function(x) gsub("\"","",gsub("gene_id ","",x[1])))
if(any(duplicated(id))) stop("non-unique gene names")
rownames(tmptab) <- id
tmptab[[ncol(tmptab)+1]] <- geneIndex
fn<-paste0(prefix,".genes.txt")
write.table(tmptab,fn,col.names=FALSE,quote=FALSE,sep="\t")
system(paste("gzip",fn))
       
###Make trascript table and get a unique ID
tmptab <- tab[txIndex,]
##get unique ID
id <- strsplit(tmptab[,9],";")
id <- sapply(id,function(x) gsub("\"","",gsub(" transcript_id ","",x[2])))
if(any(duplicated(id))) stop("non-unique transcript names")
rownames(tmptab) <- id
tmptab[[ncol(tmptab)+1]] <- txIndex
fn<-paste0(prefix,".transcripts.txt")
write.table(tmptab,fn,col.names=FALSE,quote=FALSE,sep="\t")
system(paste("gzip",fn))


###Make exon table and get a unique ID
###We do not have a way to create unique exons yet
###so no exon table yet
## tmptab <- tab[exonIndex,]
## id <- strsplit(tmptab[,9],";")





