params = RunParameters("Mayank-de-novo","menu6","$cwd/website",25,"median",0.01,"pcit",0.95,"empirical_dist",false,false,false)

#Read in raw counts (cached)
raw_counts_file = "$cwd/output/cache/$(params.test_name)_raw_counts.jld"
if (isfile(raw_counts_file))
	raw_counts = JLD.load(raw_counts_file,"raw counts")
else
	### Read in each set- generate (condensed) raw_counts and data matrix 
	raw_counts=read_count_data("data/mayank-de-novo/isoforms",method="expected_count");
	#raw_counts=RSEM.read_count_data("data/mayank-per-transcript/isoforms",method="expected_count");
	
	#filtering out into types
	code_counts=filter_count_data("data/mayank-de-novo/code-hits.list",raw_counts)
	noncode_counts=filter_count_data("data/mayank-de-novo/non-code-hits.list",raw_counts)	
	insertcols!(code_counts,"transcript_type"=>"coding")
	insertcols!(noncode_counts,"transcript_type"=>"noncoding")
	
	#merging back into one dataframe, with 
	raw_counts=outerjoin(code_counts,noncode_counts,on = intersect(names(code_counts),names(noncode_counts)))
	
	
	#boxplot(raw_counts,"raw_data_boxplot.svg")
	
	## Condense - merge polyA- and polyA+ counts for the same sample
	condensed=raw_counts[!,[1:13;26]]
	for i in 2:13
		condensed[!,i]=raw_counts[!,i]+raw_counts[!,i+12]
	end
	raw_counts=condensed
	JLD.save(raw_counts_file,"raw counts",raw_counts)
end

##run code to generate plots and figures for website 
webpage_construction(raw_counts,params)
