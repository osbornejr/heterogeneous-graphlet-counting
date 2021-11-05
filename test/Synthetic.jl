
params = RunParameters("Synthetic","menu1","$cwd/website",25,"upper_quartile",0.025,"pidc",0.95,"empirical_dist_zero",1000,true,false,true,true)


##For now we run this test with dummy GSE68559 data. The Synthetic output cache has been preprepared with the GSE runs up until the network analysis steps (obviously this could work a lot 
#Read in raw counts (cached)
raw_counts_file = "$cwd/output/cache/GSE68559_raw_counts.jld"
if (isfile(raw_counts_file))
	raw_counts = JLD.load(raw_counts_file,"raw counts")
else
	samples = CSV.read.(filter(x->occursin(".txt",x),readdir("$cwd/data/GSE68559_RAW",join=true)))
	sample_names = replace.(filter(x->occursin(".txt",x),readdir("$cwd/data/GSE68559_RAW")),"_isoforms_expr.txt"=>"").*" data"
	transcript_names = select(samples[1],1)
	gene_names = select(samples[1],4)
	gene_short_names = select(samples[1],5)
	raw_counts = rename!(hcat(transcript_names,select.(samples,Symbol("FPKM"))...,makeunique=true,gene_names,gene_short_names),["transcript_id";sample_names;"gene_id";"gene_short_id"])
	raw_counts.transcript_type = replace(x-> occursin("lnc",x) ? "noncoding" : "coding",raw_counts.transcript_id)
	JLD.save(raw_counts_file,"raw counts",raw_counts)
end

##run code to generate plots and figures for website 
webpage_construction(raw_counts,params)
