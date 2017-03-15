# Triagen
Strategy for triaging genomic variants

Attempt to automatically triage variants into those that we think are importnat, and those that we think are artefact/unimportant.

General workflow:
	1 - subset to regions of interest
	2 - remove variants that fail quality filters
	3 - filter variants believed to be artefact (germline/unidirectional)
	4 - filter variants that are common in population (ExAC)
	5 - cross reference with previous datasets (Genie,CIVIC,Sanger)
	6 - cross reference with TCGA
	7 - refine importance with predicted pathogenicity (CADD)

./R/filterNonCoding.R performs steps 1 and 2
./R/prioritiseVariants performs steps 3-7
