###read in the feature annotation file
ann=read.delim("gencode.v16.annotation.genes.gtf",skip=5,header=FALSE,stringsAsFactors=FALSE,colClasses=c("character","character","character","integer","integer","character","character","character"),comment.char="",quote="")

###read in the quantification file
y=read.delim("quantification.txt",stringsAsFactors=FALSE,row.names=1)

###read in the sample annotation file
pd=read.delim("sampleInfo.txt",stringsAsFactors=FALSE)

###get the feautre IDs from the gtf
id <- strsplit(ann[,9],";")
id <- sapply(id,function(x) gsub("\"","",gsub("feature_id ","",x[1])))


###These have to be TRUE
###sample names match?
print(identical(colnames(y),pd[,1]))
###feature names match?
print(identical(rownames(y),id))



  
