annot2vcf = function(data,outFile,chromCol=5,posCol=6,idCol=4,refCol=7,altCol=8,qualCol=9,filterCol=10,split=1)
	{
	# create vcf
	vcf = data[,c(chromCol,posCol,idCol,refCol,altCol,qualCol,filterCol)]
	newHeads = c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER")
	vcf = rbind(newHeads,vcf)
	# create vcf headers
	headers = c('##fileformat=VCFv4.0',
		paste0('##fileDate=',gsub("-","",Sys.Date())),
		'##FILTER=<ID=PASS,Description="Passed variant">')
	headFillers = matrix("",ncol=ncol(vcf)-1,nrow=length(headers))
	headers = cbind(headers,headFillers)
	colnames(headers) = colnames(vcf)
	# combined headers and data
	if(split>1)
		{
		# get file names
		ending = rev(strsplit(outFile,split="[.]")[[1]])[1]
		fileBase = gsub(paste0("[.]",ending),"",outFile)
		outFiles = paste0(fileBase,"-split-",1:split,".",ending)
		# split data
		binSize = ceiling(nrow(vcf)/split)
		cohorts = 1:nrow(vcf)
		cohorts = split(cohorts, ceiling(seq_along(cohorts)/binSize))
		vcfs = lapply(cohorts,FUN=function(x) rbind(headers,vcf[x,,drop=FALSE]))
		# write vcfs
		sapply(1:length(vcfs),FUN=function(x) write.table(vcfs[[x]],sep="\t",file=outFiles[x],row.names=FALSE,col.names=FALSE,quote=FALSE))
		} else {
		vcf = rbind(headers,vcf)
		# write vcf
		write.table(vcf,sep="\t",file=outFile,row.names=FALSE,col.names=FALSE,quote=FALSE)
		}
	}


#data = read.table("MPNST_good_Pind_Cave.txt",as.is=TRUE)
#annot2vcf(data,"MPNST_good_Pind_Cave.vcf")
