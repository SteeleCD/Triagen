#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
expectedNargs = 23
nArgs = length(args)
# test arguments: if not, return an error
if (!nArgs%in%c(1,4,expectedNargs)) 
	{
	stop("Incorrect number of arguments", call.=FALSE)
	} else {
		if(nArgs==1) 
			{
			#scriptDir=dirname(sys.frame(1)$ofile)
#scriptDir <- getSrcDirectory(function(dummy) {dummy})
#print(scriptDir)
			#setwd(scriptDir)
			args[2:4] = c("~/Dropbox/PostDoc/Rlibs/Triagen/data/knownMuts/genie_known_muts.bed",
				"~/Dropbox/PostDoc/Rlibs/Triagen/data/knownMuts/civic_variants_03022017.tsv",
				"~/Dropbox/PostDoc/Rlibs/Triagen/data/knownMuts/Sanger_drivers.csv")
			} 
		if(nArgs!=expectedNargs) 
			{
			# default column headings
			args[5:expectedNargs] = c(
				"Chrom",				# chromosome
				"Start_Position",			# start
				"End_Position",				# end
				"No.of.Occurances.in.TCGA.ICGC",	# dbCountCol 
				"ExAC_AF",				# exacCountCol 
				"Hugo_Symbol",				# geneCol 
				"HGVSc",				# variantCol 
				"IMPACT",				# impactCol 
				"CADD_PHRED",				# caddCol
				"CLIN_SIG",				# clinsigCol
				"Variant_Classification",		# variantClass
				"Sample",				# patientCol
				"Ref",					# reference
				"Alt",					# alternative
				"Type",					# variant type
				"Sub",					# substition name
				TRUE,					# bool: run unidirectional filter
				FALSE,					# bool: run germline filter
				TRUE					# bool: run CADD filter
				)
			}
		}

# running options
doUni = args[21]
doGermline = args[22]
doCADD = args[23]

classifyVariant = function(chrom, 	# chromosome
			posStart, 	# start position
			posEnd, 	# end position
			ref,		# ref seq
			alt,		# alt seq
			cDNA,		# cDNA for variant
			nPatients, 	# count of variant in patients
			tcgaCount=NULL, # count of variant in TCGA and ICGC
			exacCount, 	# count of variant in exac
			impact, 	# VEP impact
			caddScore, 	# CADD score
			genie,		# genie hotspot mutations
			clinsig,	# clinical significance
			variantClass,	# variant classification
			civic,		# civic databse
			sanger,		# sanger database
			unidirectionalFilter, # unidirectionalFilter
			germlineFilter, germlineThresh=0.015)	# germlineFilter
	{
	# civic
	matchIndex = which(paste0(civic$chromosome)==paste0(chrom)&
			civic$start==posStart&
			civic$stop==posEnd&
			paste0(civic$reference_bases)==paste0(ref)&
			paste0(civic$variant_bases)==paste0(alt))
	if(length(matchIndex)>0) return(c("High_confidence","CIVIC"))
	# sanger
	matchIndex = which(paste0(sanger$Chr)==paste0(chrom)&
			sanger$Posn==posStart&
			paste0(sanger$cDNA)==paste0(cDNA))
	if(length(matchIndex)>0) return(c("High_confidence","Sanger"))
	# genie
	variant = data.frame(chrom=chrom,start=as.numeric(posStart),end=as.numeric(posEnd))
	variant = as(variant,"GRanges")
	overlaps = findOverlaps(genie,variant)
	if(length(overlaps)>0) 
		{
		return(c("High_confidence","Genie"))
		}

	# unidirectional
	if(doUni)
		{
		if(unidirectionalFilter==0)
			{
			return(c("Unreliable","Unidirectional"))
			}
		}
	# germline filter
	if(doGermline)
		{
		if(germlineFilter>germlineThresh)
			{
			return(c("Unreliable","Germline"))
			}
		}
	# silent mutations for driver analysis
	if(variantClass=="Silent") return(c("silent","silent"))
	# exac
	if(!is.na(exacCount)) 
		{
		return(c("Unreliable","ExAC"))
		}
	# TCGA/ICGC
	if(!is.null(tcgaCount))
		{
		if(tcgaCount==0)
			{
			category = triageCurrentObs(nSamples,nPatients)
			if(doCADD) category = triagePathPred(caddScore,impact,clinsig,category)
			return(c(paste0("Novel_",category),"TCGA/CADD"))
			} else {
			category = triageCurrentObs(nSamples,nPatients)
			if(doCADD) category = triagePathPred(caddScore,impact,clinsig,category)
			return(c(paste0("Observed_",category),"TCGA/CADD"))
			}
		} else {
		category = triageCurrentObs(nSamples,nPatients)
		if(doCADD) category = triagePathPred(caddScore,impact,clinsig,category)
		return(c(category,"CADD"))
		}
	# something has gone wrong if reach here
	print("error")
	return("error")
	}

triageCurrentObs = function(nSamples,nPatients)
	{
	# more than one sample, more than one patient
	if(nPatients>1) return("low_confidence")
	# more than one sample, all in same patient
	return("ultra_low_confidence")
	}

triagePathPred = function(CADD,impact,clinsig,designation)
	{
	lowFlag = grepl("ultra",designation)	
	# increase priority of high impact 
	if(lowFlag&((CADD>30&!is.na(CADD))|impact=="HIGH"|clinsig=="pathogenic")) return(gsub("ultra_","",designation))
	# decrease priority of low impact 
	if(!lowFlag&((CADD<20&!is.na(CADD))&impact%in%c("LOW","MODIFIER"))) return(paste0("ultra_",designation))
	return(designation)
	}

# column names for columns of interest
print("set arguments")
chromCol = args[5]
posStartCol = args[6]
posEndCol = args[7]
tcgaCountCol = args[8]
exacCountCol = args[9]
geneCol = args[10]
variantCol = args[11] 
impactCol = args[12]
caddCol = args[13]
clinsigCol = args[14]
variantClassCol = args[15]
patientCol = args[16]
refCol = args[17]
altCol = args[18]
variantTypeCol = args[19]
# variable names
subName = args[20]
# data files
dataFile = args[1]
genieFile = args[2]
civicFile = args[3]
sangerFile = args[4]
# load data
print("read data")
data = read.table(dataFile,head=TRUE,sep="\t",comment.char="@",quote="",as.is=TRUE)
# load genie
print("read genie")
genie = read.csv(genieFile,head=FALSE,skip=10)
# convert genie to 1-based
genie[,2] = genie[,2]+1
colnames(genie)[1:3] = c("chrom","start","end")
library(GenomicFeatures)
genie = as(genie,"GRanges")
# load civic
print("read civic")
civic = read.table(civicFile,head=TRUE,sep="\t",comment.char="@",quote="")
civic = civic[which(civic$reference_bases!=""),]
# load sanger
print("read sanger")
sanger = read.csv(sangerFile,head=TRUE)
# get unidirectional filter
print("set unidirectional filter")
getMinStrand = function(data)
	{
	if(data[,variantTypeCol]==subName)
		{
		# substitutions - minimum of alternative
		alt = data[,altCol]
		return(min(data[,paste0("F",alt,"Z.Tum")],data[,paste0("R",alt,"Z.Tum")]))
		} else {
		# indels - minimum of unique calls (Pindel and BWA)
		return(min(c(data$PU.Tum,data$NU.Tum)))
		}
	}
if(doUni)
	{
	unidirectionalFlag = sapply(1:nrow(data),FUN=function(x) getMinStrand(data[x,,drop=FALSE]))
	} else {
	unidirectionalFlag = rep(NA,nrow(data))
	}

data = cbind(data,unidirectionalFlag)
uniCol = "unidirectionalFlag"
# get germline filter
print("set germline filter")
getGermline = function(data,infoMut=c("PU.Norm","NU.Norm"),infoAll=c("PR.Norm","NR.Norm"))
	{
	if(data[,variantTypeCol]==subName)
		{
		# substitutions - normal alt / normal alt & ref
		alt = data[,altCol]
	        ref = data[,refCol]
		altN = sum(as.numeric(data[,paste0(c("F","R"),alt,"Z.Norm")]))
		refN = sum(as.numeric(data[,paste0(c("F","R"),ref,"Z.Norm")]))
		return(altN/(altN+refN))
		} else {
		# indels - normal mutant / normal total
		return(sum(as.numeric(unlist(data[,infoMut])))/sum(as.numeric(unlist(data[,infoAll]))))
		}
	}
if(doGermline)
	{
	germlineFlag = sapply(1:nrow(data),FUN=function(x) getGermline(data[x,,drop=FALSE]))
	} else {
	germlineFlag = rep(NA,times=nrow(data))
	}
data = cbind(data,germlineFlag)
germCol = "germlineFlag"







# function to count the number of patients with the variant before classifying
prioritiseVariant = function(chrom,posStart,posEnd,
			gene,variant,
			tcgaCount=NULL,
			exacCount,
			impact,
			cadd,
			clinsig,
			variantClass,
			ref,
			alt,
			unidirectional,
			germline,
			patient)
	{
INDEX <<- INDEX+1
	print(paste0(INDEX,". ",gene,":",variant,",",chrom,":",posStart,"-",posEnd))
	# subset to this gene and variant
	varIndex = which(data[,geneCol]==gene&data[,variantCol]==variant)
	subData = data[varIndex,]
	# count number of patients with this variant
	patientCount = length(unique(patient[varIndex]))
	# classify variant priority
	priority = classifyVariant(
		  chrom=chrom,
		  posStart=posStart,
		  posEnd=posEnd,
		  cDNA=variant,
		  nPatients=patientCount,
                  tcgaCount=tcgaCount,
                  exacCount=exacCount,
		  impact=impact,
		  caddScore=cadd,
		  genie=genie,
		  clinsig=clinsig,
		  variantClass=variantClass,
		  sanger=sanger,
		  civic=civic,
		  ref=ref,
		  alt=alt,
		  unidirectionalFilter=unidirectional,
		  germlineFilter=germline)
	return(c(priority,patientCount))
	}

# set patient
print("set patient if missing")
if(!patientCol%in%colnames(data))
	{
	patient = rep("patient",nrow(data))
	} else {
	patient = data[,patientCol]
	}

toRun = 1:nrow(data)
# perform classification
print("run prioritisation")
INDEX<<-1
Priority = mapply(FUN=prioritiseVariant,
		    chrom=data[toRun,chromCol],
		    posStart=data[toRun,posStartCol],
		    posEnd=data[toRun,posEndCol],
                    gene=data[toRun,geneCol],
                    variant=data[toRun,variantCol],
                    tcgaCount=data[toRun,tcgaCountCol],
                    exacCount=data[toRun,exacCountCol],
		    impact=data[toRun,impactCol],
		    cadd=data[toRun,caddCol],
		    clinsig=data[toRun,clinsigCol],
		    variantClass=data[toRun,variantClassCol],
		    ref=data[toRun,refCol],
		    alt=data[toRun,altCol],
		    unidirectional=data[toRun,uniCol],
		    germline=data[toRun,germCol],
		    patient=patient
                    )
rownames(Priority) = c("priority","reason","patientCount")
priorityTable = table(Priority["priority",])
reasonTable = table(Priority["reason",])

# combine data and results
data = cbind(data,t(Priority))

# get out file name
print("set file name")
fileEnding = rev(strsplit(dataFile,split="[.]")[[1]])[1]
outFileSave = rev(strsplit(dataFile,"/")[[1]])[1]
outFile = gsub(outFileSave,"",dataFile)
# increment version
chars = strsplit(outFileSave,split="")[[1]]
if(chars[1]=="v") 
	{
	chars[2] = as.numeric(chars[2])+1
	outFileEnd = paste0(chars,collapse="") 
	} else {
	outFileEnd = outFileSave
	}
# write results
print("write results")
write.table(data,
	file=paste0(outFile,outFileEnd,"-withPriority.txt"),
	quote=FALSE,row.names=FALSE,sep="\t")














