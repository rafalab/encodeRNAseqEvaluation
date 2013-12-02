VERSION<-2
SKIP<-5 ##number of lines to skip to get the information
##get the number of rows to read in
filename<- "gencode.v16.annotation.gtf"
prefix<-"gencode.v16.annotation"
nrow=system(paste("wc -l",filename),intern=TRUE)
nrow=as.numeric(strsplit(nrow," ")[[1]][2])-SKIP

##read in table
if(!exists("tab")) tab <-  read.delim(filename,skip=SKIP,header=FALSE,stringsAsFactors=FALSE,colClasses=c("character","character","character","integer","integer","character","character","character"),comment.char="",nrow=nrow,quote="")

##split by gene, transcript, exon
Indexes <- split(1:nrow(tab),tab[,3])
geneIndex <- Indexes[["gene"]]
txIndex <- Indexes[["transcript"]]
exonIndex <- Indexes[["exon"]]

###Make gene table and get a unique ID
tmptab <- tab[geneIndex,]
###get a unique ID
id <- strsplit(tmptab[,9],";")
id <- sapply(id,function(x) gsub("\"","",gsub("gene_id ","",x[1])))
if(any(duplicated(id))) stop("non-unique gene names")
##add ninth column
tmp1<-paste0("feature_id \"",id,"\"")
tmp2<-paste("original_row",geneIndex)
tmp3<-paste(tmp1,tmp2,tmptab[,9],sep="; ")
tmptab[,9]<-tmp3


fn<-paste0(prefix,".genes.gtf")
cat(paste("##description: subset of",filename,"including only genes and adds ENCODE ERCC spikeins, version",VERSION,"\n"),file=fn)
cat("##provider: ENCODE DAC\n",file=fn,append=TRUE)
cat("##contact: rafa@jimmy.harvard.edu\n",file=fn,append=TRUE)
cat("##format: gtf\n",file=fn,append=TRUE)
cat(paste("##date:",date(),"\n"),file=fn,append=TRUE)

write.table(tmptab,fn,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
##Add spikein
spikein<-read.delim("spikein_gene.gtf",header=FALSE,stringsAsFactors=FALSE,colClasses=c("character","character","character","integer","integer","character","character","character"),comment.char="",quote="")
spikein[[2]] <- "spike_in"
write.table(spikein,fn,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
system(paste("gzip -f",fn))
id=c(id,spikein[,1])
write.table(id,file="genes.id.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)


###Make trascript table and get a unique ID
tmptab <- tab[txIndex,]
##get unique ID
id <- strsplit(tmptab[,9],";")
id <- sapply(id,function(x) gsub("\"","",gsub(" transcript_id ","",x[2])))
if(any(duplicated(id))) stop("non-unique transcript names")
##add column
tmp1<-paste0("feature_id \"",id,"\"")
tmp2<-paste("original_row",geneIndex)
tmp3<-paste(tmp1,tmp2,tmptab[,9],sep="; ")
tmptab[,9]<-tmp3

fn<-paste0(prefix,".transcripts.gtf")
cat(paste("##description: subset of",filename,"including only transcripts and adds ENCODE ERCC spikeins, version",VERSION,"\n"),file=fn)
cat("##provider: ENCODE DAC\n",file=fn,append=TRUE)
cat("##contact: rafa@jimmy.harvard.edu\n",file=fn,append=TRUE)
cat("##format: gtf\n",file=fn,append=TRUE)
cat(paste("##date:",date(),"\n"),file=fn,append=TRUE)
write.table(tmptab,fn,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
##Add spikein
spikein<-read.delim("spikein_transcript.gtf",header=FALSE,stringsAsFactors=FALSE,colClasses=c("character","character","character","integer","integer","character","character","character"),comment.char="",quote="")
spikein[[2]] <- "spike_in"
write.table(spikein,fn,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
system(paste("gzip -f",fn))
id=c(id,spikein[,1])
write.table(id,file="transcript.id.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

###Make exon table and get a unique ID
###We do not have a way to create unique exons yet
###so no exon table yet
tmptab <- tab[exonIndex,]
id <- strsplit(tmptab[,9],";")
##using Julien Lagarde suggestions
id1 <- sapply(id,function(x) gsub("\"","",gsub(" transcript_id ","",x[2])))
id2 <- sapply(id,function(x) gsub("\"","",gsub(" exon_number ","",x[9])))
id <- paste(id1,id2,sep="_")
if(any(duplicated(id))) stop("non-unique transcript names")
write.table(id,file="exon.id.txt",row.names=FALSE,col.names=FALSE)
###add column
tmp1<-paste0("feature_id \"",id,"\"")
tmp2<-paste("original_row",exonIndex)
tmp3<-paste(tmp1,tmp2,tmptab[,9],sep="; ")
tmptab[,9]<-tmp3

fn<-paste0(prefix,".exon.gtf")
cat(paste("##description: subset of",filename,"including only exons and adds ENCODE ERCC spikeins, version",VERSION,"\n"),file=fn)
cat("##provider: ENCODE DAC\n",file=fn,append=TRUE)
cat("##contact: rafa@jimmy.harvard.edu\n",file=fn,append=TRUE)
cat("##format: gtf\n",file=fn,append=TRUE)
cat(paste("##date:",date(),"\n"),file=fn,append=TRUE)
write.table(tmptab,fn,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE) 
##Add spikein
spikein<-read.delim("spikein_exon.gtf",header=FALSE,stringsAsFactors=FALSE,colClasses=c("character","character","character","integer","integer","character","character","character"),comment.char="",quote="")
spikein[[2]] <- "spike_in"
write.table(spikein,fn,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
system(paste("gzip -f",fn))
id=c(id,spikein[,1])
write.table(id,file="exon.id.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)




