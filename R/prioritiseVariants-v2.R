#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
expectedNargs = 27
nArgs = length(args)
# test arguments: if not, return an error
if (!nArgs%in%c(1,5,expectedNargs)) 
	{
	stop("Incorrect number of arguments", call.=FALSE)
	} else {
		if(nArgs==1) 
			{
			#scriptDir=dirname(sys.frame(1)$ofile)
#scriptDir <- getSrcDirectory(function(dummy) {dummy})
#print(scriptDir)
			#setwd(scriptDir)
			args[2:5] = c("~/Dropbox/PostDoc/Rlibs/Triagen/data/knownMuts/genie_mutations_2.0.0_2017-11-27.txt",
				"~/Dropbox/PostDoc/Rlibs/Triagen/data/knownMuts/civic_variants_2017-11-01.txt",
				"~/Dropbox/PostDoc/Rlibs/Triagen/data/knownMuts/Sanger_drivers.csv",
				"~/Dropbox/PostDoc/Rlibs/Triagen/data/knownMuts/mskcc_hotspots_2017-11-27.txt")
			} 
		if(nArgs!=expectedNargs) 
			{
			# default column headings
			args[6:expectedNargs] = c(
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
				TRUE,					# bool: run CADD filter
				"sanger",				# germline method
				"n_ref_count",				# normal reference count
				"n_alt_count"				# normal alternative count
				)
			}
		}



# function to prioritise variants based on a number of factors
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
			germlineFilter, germlineThresh=0.015, # germlineFilter
			doUni=TRUE,doGermline=FALSE,doCADD=TRUE, # running options
			exacThresh=0.0004) # exac threshold
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
		print(germlineFilter)
		if(!is.na(germlineFilter))
			{
			if(germlineFilter>germlineThresh)
				{
				return(c("Unreliable","Germline"))
				}
			}
		}
	# silent mutations for driver analysis
	if(variantClass=="Silent") return(c("silent","silent"))
	# exac
	if(!is.na(exacCount)&exacCount>exacThresh) 
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

# function to triage base on numbers of observations in samples/patients
triageCurrentObs = function(nSamples,nPatients)
	{
	# more than one sample, more than one patient
	if(nPatients>1) return("low_confidence")
	# more than one sample, all in same patient
	return("ultra_low_confidence")
	}

# function to triage based on predicted pathogenicity
triagePathPred = function(CADD,impact,clinsig,designation)
	{
	lowFlag = grepl("ultra",designation)	
	# increase priority of high impact 
	if(lowFlag&((CADD>30&!is.na(CADD))|impact=="HIGH"|(!is.na(clinsig)&clinsig=="pathogenic")))
		{
		return(gsub("ultra_","",designation))
		}
	# decrease priority of low impact 
	if(!lowFlag&((CADD<20&!is.na(CADD))&impact%in%c("LOW","MODIFIER")))
		{
		return(paste0("ultra_",designation))
		}
	return(designation)
	}


# function to check validity of columns
checkCol = function(data,column)
	{
	check = column%in%colnames(data)
	if(!check) stop(paste0("Column ",column," missing from data"))
	if(all(is.na(data[,column]))) print(paste0("Warning: all data values are NA in column ",column))
	}

# function to check files
checkFile = function(fileName)
	{
	# check if file exists
	if(!file.exists(fileName)) stop(paste0("File ",fileName," does not exist"))
	# check if file is empty
	info = file.info(fileName)
	if(info$size==0) stop(paste0("File ",fileName," is empty"))
	}

# get out filename
getOutFile = function(filename)
	{
	# get out file name
	print("set file name")
	fileEnding = rev(strsplit(filename,split="[.]")[[1]])[1]
	outFileSave = rev(strsplit(filename,"/")[[1]])[1]
	outFile = gsub(outFileSave,"",filename)
	# increment version
	chars = strsplit(outFileSave,split="")[[1]]
	if(chars[1]=="v") 
		{
		chars[2] = as.numeric(chars[2])+1
		outFileEnd = paste0(chars,collapse="") 
		} else {
		outFileEnd = outFileSave
		}
	# return filename
	return(paste0(outFile,outFileEnd,"-withPriority.txt"))
	}

# function to get germline filter
getGermlineSanger = function(data,infoMut=c("PU.Norm","NU.Norm"),infoAll=c("PR.Norm","NR.Norm"),
		variantTypeCol,altCol,refCol,subName)
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


# function to get unidirectional filter
getMinStrand = function(data,variantTypeCol,altCol,subName)
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

# funciton to set up prioritise
setupPrioritise = function(dataFile,genieFile,civicFile,sangerFile,chromCol,
		posStartCol,posEndCol,tcgaCountCol,exacCountCol,geneCol,
		variantCol,impactCol,caddCol,clinsigCol,variantClassCol,
		patientCol,refCol,altCol,variantTypeCol,subName,
		doUni=TRUE,doGermline=FALSE,doCADD=TRUE,
		germlineMeth="sanger",germlineRefCountCol=NA,germlineAltCountCol=NA)
	{	
	# check input files
	sapply(c(dataFile,
		genieFile,
		civicFile,
		sangerFile),FUN=checkFile)
	# load data
	print("read data")
	data = read.table(dataFile,head=TRUE,sep="\t",comment.char="@",quote="",as.is=TRUE)
	# check data columns
	print("check data")
	sapply(c(chromCol,
		posStartCol,
		posEndCol,
		tcgaCountCol,
		exacCountCol,
		geneCol,
		variantCol,
		impactCol,
		clinsigCol,
		variantClassCol,
		refCol,
		altCol,
		variantTypeCol),FUN=function(x) checkCol(data,x))
	if(doCADD) checkCol(data,caddCol) 
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
	if(doUni)
		{
		unidirectionalFlag = sapply(1:nrow(data),FUN=function(x)
			{
			getMinStrand(data[x,,drop=FALSE],
				variantTypeCol=variantTypeCol,
				altCol=altCol,subName=subName)
			})
		} else {
		unidirectionalFlag = rep(NA,nrow(data))
		}
	data = cbind(data,unidirectionalFlag)
	uniCol = "unidirectionalFlag"
	# get germline filter
	print("set germline filter")
	if(doGermline)
		{
		if(germlineMeth=="sanger")
			{
			germlineFlag = sapply(1:nrow(data),
				FUN=function(x) getGermlineSanger(data[x,,drop=FALSE],
						variantTypeCol=variantTypeCol,
						altCol=altCol,
						refCol=refCol,
						subName=subName))
			} else {
			if(!is.na(germlineRefCountCol)&!is.na(germlineAltCountCol))
				{
				denom = (as.numeric(data[,germlineAltCountCol])+as.numeric(data[,germlineRefCountCol]))
				germlineFlag = as.numeric(data[,germlineAltCountCol])/denom
				} else {
				germlineFlag = rep(NA,times=nrow(data))
				}
			}
		} else {
		germlineFlag = rep(NA,times=nrow(data))
		}
	data = cbind(data,germlineFlag)
	germCol = "germlineFlag"
	# set patient
	print("set patient if missing")
	if(!patientCol%in%colnames(data))
		{
		patient = rep("patient",nrow(data))
		} else {
		patient = data[,patientCol]
		}
	# return data
	return(list(data=data,civic=civic,sanger=sanger,
		genie=genie,uniCol=uniCol,germCol=germCol,
		patient=patient))
	}


# function to count the number of patients with the variant before prioritising
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
			patient,
			geneCol,
			variantCol,
			civic,
			genie,
			sanger,
			data,
			doUni,doGermline,doCADD
			)
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
		  germlineFilter=germline,
		  doUni=doUni,doGermline=doGermline,doCADD=doCADD)
	return(c(priority,patientCount))
	}



# function to run priorities in full (pipeline)
runPriorities  = function(dataFile,genieFile,civicFile,sangerFile,
		chromCol,posStartCol,posEndCol,tcgaCountCol,exacCountCol,
		geneCol,variantCol,impactCol,caddCol,clinsigCol,variantClassCol,
		patientCol,refCol,altCol,variantTypeCol,subName,
		doUni=TRUE,doGermline=FALSE,doCADD=TRUE,
		germlineMeth="sanger",germlineRefCountCol=NA,germlineAltCountCol=NA,
		toRun="all",doWrite=TRUE)
	{
	# read and check data
	data = setupPrioritise(dataFile,genieFile,civicFile,sangerFile,
		chromCol,posStartCol,posEndCol,tcgaCountCol,exacCountCol,
		geneCol,variantCol,impactCol,caddCol,clinsigCol,variantClassCol,
		patientCol,refCol,altCol,variantTypeCol,subName,
		doUni,doGermline,doCADD,germlineMeth,germlineRefCountCol,germlineAltCountCol)
	DATA <<- data
	# which variants to run
	if(toRun=="all") toRun = 1:nrow(data$data)
	# perform classification
	print("run prioritisation")
	INDEX<<-1
	Priority = mapply(FUN=prioritiseVariant,
		    chrom=data$data[toRun,chromCol],
		    posStart=data$data[toRun,posStartCol],
		    posEnd=data$data[toRun,posEndCol],
                    gene=data$data[toRun,geneCol],
                    variant=data$data[toRun,variantCol],
                    tcgaCount=data$data[toRun,tcgaCountCol],
                    exacCount=data$data[toRun,exacCountCol],
		    impact=data$data[toRun,impactCol],
		    cadd=data$data[toRun,caddCol],
		    clinsig=data$data[toRun,clinsigCol],
		    variantClass=data$data[toRun,variantClassCol],
		    ref=data$data[toRun,refCol],
		    alt=data$data[toRun,altCol],
		    unidirectional=data$data[toRun,data$uniCol],
		    germline=data$data[toRun,data$germCol],
		    MoreArgs=list(geneCol=geneCol,variantCol=variantCol,
				civic=data$civic,sanger=data$sanger,
				genie=data$genie,data=data$data,
				patient=data$patient,
				doUni=doUni,doGermline=doGermline,doCADD=doCADD)
                    )
	rownames(Priority) = c("priority","reason","patientCount")
	priorityTable = table(Priority["priority",])
	reasonTable = table(Priority["reason",])
	# combine data and results
	data = cbind(data$data,t(Priority))
	# write results
	if(doWrite)
		{
		outFile = getOutFile(dataFile)
		print("write results")
		write.table(data,
			file=outFile,
			quote=FALSE,row.names=FALSE,sep="\t")
		}
	return(data)
	}

# collate arguments into a list
collateArgs = function(args)
	{
	print("set arguments")
	args = list(
		# data files
		dataFile = args[1],
		genieFile = args[2],
		civicFile = args[3],
		sangerFile = args[4],
		mskccFile = args[5],
		# column names
		chromCol = args[6],
		posStartCol = args[7],
		posEndCol = args[8],
		tcgaCountCol = args[9],
		exacCountCol = args[10],
		geneCol = args[11],
		variantCol = args[12], 
		impactCol = args[13],
		caddCol = args[14],
		clinsigCol = args[15],
		variantClassCol = args[16],
		patientCol = args[17],
		refCol = args[18],
		altCol = args[19],
		variantTypeCol = args[20],
		# variable names
		subName = args[21],
		# running options
		doUni = args[22],
		doGermline = args[23],
		doCADD = args[24],
		germlineMeth = args[25],
		germlineRefCountCol = args[26],
		germlineAltCountCol = args[27]
		)
	return(args)
	}



# prioritise
prioritiseVars = function(args)
	{
	argsList = collateArgs(args)
	do.call(runPriorities,argsList)
	}




# function to run priorities in full (pipeline), but only once per variant
runPrioritiesOnce  = function(dataFile,genieFile,civicFile,sangerFile,
		chromCol,posStartCol,posEndCol,tcgaCountCol,exacCountCol,
		geneCol,variantCol,impactCol,caddCol,clinsigCol,variantClassCol,
		patientCol,refCol,altCol,variantTypeCol,subName,
		doUni=TRUE,doGermline=FALSE,doCADD=TRUE,
		germlineMeth="sanger",germlineRefCountCol=NA,germlineAltCountCol=NA,
		toRun="all",doWrite=TRUE,outFile=NULL)
	{
	# read and check data
	data = setupPrioritise(dataFile,genieFile,civicFile,sangerFile,
		chromCol,posStartCol,posEndCol,tcgaCountCol,exacCountCol,
		geneCol,variantCol,impactCol,caddCol,clinsigCol,variantClassCol,
		patientCol,refCol,altCol,variantTypeCol,subName,
		doUni,doGermline,doCADD,germlineMeth,germlineRefCountCol,germlineAltCountCol)
	DATA <<- data
	# which variants to run
	if(toRun=="all") toRun = 1:nrow(data$data)
	# setup so only run once for each variant
	dictionary = paste0(data$data[toRun,chromCol],":",
			data$data[toRun,posStartCol],":",
			data$data[toRun,posEndCol],":",
			data$data[toRun,refCol],":",
			data$data[toRun,altCol],":",
			data$data[toRun,geneCol],":",
			data$data[toRun,variantCol]
			)
	dictUnique = unique(dictionary)
	dictIndices = sapply(dictUnique,FUN=function(x) which(dictionary==x))
	dictStarts = sapply(dictIndices,FUN=function(x) x[1])
	rm(dictionary)
	rm(dictUnique)
	gc()
	# perform classification once per variant
	print("run prioritisation")
	INDEX<<-1
	Priority = mapply(FUN=prioritiseVariant,
		    chrom=data$data[toRun[dictStarts],chromCol],
		    posStart=data$data[toRun[dictStarts],posStartCol],
		    posEnd=data$data[toRun[dictStarts],posEndCol],
                    gene=data$data[toRun[dictStarts],geneCol],
                    variant=data$data[toRun[dictStarts],variantCol],
                    tcgaCount=data$data[toRun[dictStarts],tcgaCountCol],
                    exacCount=data$data[toRun[dictStarts],exacCountCol],
		    impact=data$data[toRun[dictStarts],impactCol],
		    cadd=data$data[toRun[dictStarts],caddCol],
		    clinsig=data$data[toRun[dictStarts],clinsigCol],
		    variantClass=data$data[toRun[dictStarts],variantClassCol],
		    ref=data$data[toRun[dictStarts],refCol],
		    alt=data$data[toRun[dictStarts],altCol],
		    unidirectional=data$data[toRun[dictStarts],data$uniCol],
		    germline=data$data[toRun[dictStarts],data$germCol],
		    MoreArgs=list(geneCol=geneCol,variantCol=variantCol,
				civic=data$civic,sanger=data$sanger,
				genie=data$genie,data=data$data,
		    		patient=data$patient,
				doUni=doUni,doGermline=doGermline,doCADD=doCADD)
                    )
	rownames(Priority) = c("priority","reason","patientCount")
	# fill out object with duplicate variants
	out = matrix(NA,ncol=3,nrow=length(toRun))
	colnames(out) = c("priority","reason","patientCount")
	for(i in 1:length(dictIndices))
		{
		for(j in 1:length(dictIndices[[i]]))
			{
			out[dictIndices[[i]][j],] = Priority[,i]
			}
		}
	# tables
	priorityTable = table(out[,"priority"])
	reasonTable = table(out[,"reason"])
	# combine data and results
	data = cbind(data$data,out)
	# write results
	if(doWrite)
		{
		if(is.null(outFile)) outFile = getOutFile(dataFile)
		print("write results")
		write.table(data,
			file=outFile,
			quote=FALSE,row.names=FALSE,sep="\t",na="")
		}
	return(data)
	}

# prioritise
prioritiseVarsOnce = function(args)
	{
	argsList = collateArgs(args)
	do.call(runPrioritiesOnce,argsList)
	}

# actually run prioritisation
prioritiseVarsOnce(args)


