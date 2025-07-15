using DataFrames,JLD2,ProjectFunctions,CSV,DataPreprocessing
pwd = ENV["PWD"]
### Read in all counts
input_dir ="$pwd/data/mayank-de-novo"
@info "Pulling in all raw counts..."
raw_counts=DataPreprocessing.read_count_data("$input_dir/isoforms",method="expected_count");
@info "found $(size(raw_counts)[1]) transcripts at $input_dir..."   
#raw_counts=RSEM.read_count_data("data/mayank-per-transcript/isoforms",method="expected_count");

#filtering out into types
@info "filtering out to just classified transcripts..."
code_counts=DataPreprocessing.filter_count_data("$pwd/data/mayank-de-novo/code-hits.list",raw_counts)

@info "found $(size(code_counts)[1]) coding transcripts..."   
noncode_counts=DataPreprocessing.filter_count_data("$pwd/data/mayank-de-novo/non-code-hits.list",raw_counts)	

@info "found $(size(noncode_counts)[1]) coding transcripts..."   

insertcols!(code_counts,"transcript_type"=>"coding")
insertcols!(noncode_counts,"transcript_type"=>"noncoding")
#merging back into one dataframe, with 

@info "removing non-classified transcripts..."   
raw_counts=outerjoin(code_counts,noncode_counts,on = intersect(names(code_counts),names(noncode_counts)))

@info "merging polyA+ and polyA- counts..."

condensed=raw_counts[!,[1:13;26]]
for i in 2:13
	condensed[!,i]=raw_counts[!,i]+raw_counts[!,i+12]
end 
raw_counts=condensed
rename!(raw_counts,replace.(names(raw_counts),"polyA+_" => ""))

file = "$pwd/output/cache/mayank-merged/raw_counts.jld2"

@info "Saving raw counts to $file"
cache_save(file,"raw counts"=>raw_counts)
