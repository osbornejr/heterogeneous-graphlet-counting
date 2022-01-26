using RCall,DataFrames,JLD2 
#params = ProjectFunctions.RunParameters("GSE68559","menu1","$cwd/website",25,"upper_quartile",0.025,"pidc",0.95,"empirical_dist_zero",1000,true,false,true,true)

#Read in raw counts (cached)
raw_counts_file = "$cwd/output/cache/$(params["test_name"])_raw_counts.jld2"
if (isfile(raw_counts_file))
    raw_counts = cache_load(raw_counts_file,"raw counts")
else
	samples = CSV.read.(filter(x->occursin(".txt",x),readdir("$cwd/data/GSE68559_RAW",join=true)))
	sample_names = replace.(filter(x->occursin(".txt",x),readdir("$cwd/data/GSE68559_RAW")),"_isoforms_expr.txt"=>"").*" data"
	transcript_names = select(samples[1],1)
	gene_names = select(samples[1],4)
	gene_short_names = select(samples[1],5)
	raw_counts = rename!(hcat(transcript_names,select.(samples,Symbol("FPKM"))...,makeunique=true,gene_names,gene_short_names),["transcript_id";sample_names;"gene_id";"gene_short_id"])
	raw_counts.transcript_type = replace(x-> occursin("lnc",x) ? "noncoding" : "coding",raw_counts.transcript_id)
	cache_save(raw_counts_file,"raw counts"=>raw_counts)
end

## find transcript types via biomaRt (NOTE this should be inside the raw counts cache section above, but currently input data has vanished from NeCTAR! It is all cached here thankfully, so for now we use this workaround.

biomart_modification = true


if(biomart_modification = true)
    biomart_raw_counts_file = "$cwd/output/cache/$(params["test_name"])_raw_counts_biomart.jld2"
	if (isfile(biomart_raw_counts_file))
        raw_counts = cache_load(biomart_raw_counts_file,"raw counts")
	else
		#restart R session	
		R"""
		sapply(names(sessionInfo()$otherPkgs),function(pkg) detach(paste0('package:',pkg),character.only =T,force = T));rm(list=ls())
		"""
		transcripts = raw_counts.transcript_id
		@rput transcripts
		R"""
		library(biomaRt)
		library(httr)
		library(edgeR)
		library(tidyverse)
		transcripts_trimmed = sapply(transcripts,tools::file_path_sans_ext)
		## connect to biomart
		set_config(config(ssl_verifypeer = 0L))
		ensembl_version = "current"	
		if (ensembl_version=="current")
			{
			##mirrors to try: "useast" "uswest" "asia"
			ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",mirror="uswest", dataset = "hsapiens_gene_ensembl") 
			} else
			{
			ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",version=ensembl_version) 
			}
		
		transcripts_biotypes = getBM(attributes=c("ensembl_transcript_id","transcript_biotype"),"ensembl_transcript_id",transcripts_trimmed, mart = ensembl,useCache=FALSE)
		transcript_types= transcripts_biotypes[match(transcripts_trimmed,transcripts_biotypes[[1]]),2]
		"""
		
		@rget transcript_types
		## get rid of missing values and add to raw counts data
		raw_counts.biomart_type = replace(transcript_types,missing=>"none")
		##cut down raw counts to only those with coding or lncRNA labels (as well as those explicitly labelled lncRNA in input data).
		raw_counts = vcat(filter(:transcript_type=>x->x=="noncoding",raw_counts),filter(:biomart_type=>x->x in ["lncRNA","protein_coding"],raw_counts))
		## set transcripts that biomart identifies as lncRNA to noncoding
		raw_counts.transcript_type = @. ifelse(raw_counts.biomart_type=="lncRNA","noncoding",raw_counts.transcript_type)
		cache_save(biomart_raw_counts_file,"raw counts"=>raw_counts)
	end
end
##run code to generate plots and figures for website 
#webpage_construction(raw_counts,params)
