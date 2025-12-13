using RCall,DataFrames,JLD2,ProjectFunctions,CSV 
#params = ProjectFunctions.RunParameters("GSE68559","menu1","$cwd/website",25,"upper_quartile",0.025,"pidc",0.95,"empirical_dist_zero",1000,true,false,true,true)

#### This is a dataset specific function that will generate the raw counts file associated with GSE68559

cwd = ENV["PWD"]
### get data from website if necessary
if !(isdir("$cwd/data/GSE68559/GSE68559_RAW"))
    run(`mkdir -p $cwd/data/GSE68559`)
    run(`wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE68nnn/GSE68559/suppl/GSE68559_RAW.tar`)
    run(`tar -xf GSE68559_RAW.tar --one-top-level`)
    run(`mv GSE68559_RAW $cwd/data/GSE68559/`)
    for file in readdir("$cwd/data/GSE68559/GSE68559_RAW")
        if occursin("het.vcf",file)
            run(`rm $cwd/data/GSE68559/GSE68559_RAW/$file`)
        else
            run(`gunzip $cwd/data/GSE68559/GSE68559_RAW/$file`)
        end
    end
    run(`rm GSE68559_RAW.tar`)
end



#Read in raw counts 
@info "Reading in isoform data"
samples = CSV.read.(filter(x->occursin(".txt",x),readdir("$cwd/data/GSE68559/GSE68559_RAW",join=true)),Ref(DataFrame))
sample_names = replace.(filter(x->occursin(".txt",x),readdir("$cwd/data/GSE68559/GSE68559_RAW")),"_isoforms_expr.txt"=>"").*" data"
transcript_names = select(samples[1],1)
gene_names = select(samples[1],4)
gene_short_names = select(samples[1],5)
raw_counts = rename!(hcat(transcript_names,select.(samples,Symbol("FPKM"))...,makeunique=true,gene_names,gene_short_names),["transcript_id";sample_names;"gene_id";"gene_short_id"])
##pre add this column, denoting just those that are explicitly lnc labelled as noncoding
raw_counts.transcript_type = replace(x-> occursin("lnc",x) ? "noncoding" : "coding",raw_counts.transcript_id)
###2025 add: rename duplicates here so that they match R make.unique
##for simplicity we will take them directly from the biotypes file that is now generated in an R script
transcript_biotypes = CSV.read("data/GSE68559/transcript_biotypes.txt",DataFrame)
raw_counts.transcript_id = transcript_biotypes.TranscriptID

## find transcript types via biomaRt (NOTE this should be inside the raw counts cache section above, but currently input data has vanished from NeCTAR! It is all cached here thankfully, so for now we use this workaround.
#@info "Checking transcript types via biomaRt"
#restart R session	
#R"""
#sapply(names(sessionInfo()$otherPkgs),function(pkg) detach(paste0('package:',pkg),character.only =T,force = T));rm(list=ls())
#"""
#transcripts = raw_counts.transcript_id
#@rput transcripts
#R"""
#library(biomaRt,quietly = T)
#library(httr,quietly = T)
#library(edgeR,quietly = T)
#library(tidyverse,quietly = T)
#transcripts_trimmed = sapply(transcripts,tools::file_path_sans_ext)
### connect to biomart
#set_config(config(ssl_verifypeer = 0L))
#ensembl_version = "109"	
#if (ensembl_version=="current")
#{
#    ##mirrors to try: "useast" "uswest" "asia"
#    ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",mirror="uswest", dataset = "hsapiens_gene_ensembl") 
#} else
#{
#    ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",version=ensembl_version) 
#}
#
#transcripts_biotypes = getBM(attributes=c("ensembl_transcript_id","transcript_biotype"),"ensembl_transcript_id",transcripts_trimmed, mart = ensembl,useCache=FALSE)
#transcript_types= transcripts_biotypes[match(transcripts_trimmed,transcripts_biotypes[[1]]),2]
#"""
#
#@rget transcript_types
## get rid of missing values and add to raw counts data
##2025 add: now we get the biomart types externally, loaded viave transcrpt biotypes file
#raw_counts.biomart_type = replace(transcript_types,missing=>"none")
raw_counts.biomart_type = transcript_biotypes
##cut down raw counts to only those with coding or lncRNA labels (as well as those explicitly labelled lncRNA in input data).
raw_counts = vcat(filter(:transcript_type=>x->x=="noncoding",raw_counts),filter(:biomart_type=>x->x in ["lncRNA","protein_coding"],raw_counts))
## set transcripts that biomart identifies as lncRNA to noncoding
raw_counts.transcript_type = @. ifelse(raw_counts.biomart_type=="lncRNA","noncoding",raw_counts.transcript_type)

##fix string types to all be standard strings
for i in 1:size(raw_counts)[2]
    if (eltype(raw_counts[:,i])<:AbstractString)
        raw_counts[!,i] = String.(raw_counts[:,i])
    end
end

file = "$cwd/output/cache/GSE68559/raw_counts.jld2"

@info "Saving raw counts to $file"
cache_save(file,"raw counts"=>raw_counts)
##run code to generate plots and figures for website 
#webpage_construction(raw_counts,params)
