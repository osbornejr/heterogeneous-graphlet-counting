using Pkg
Pkg.activate(".")
using GLMakie,StatsBase
using ProjectFunctions
raw_counts,round_counts,vst_counts,clean_counts,norm_counts,processed_counts = get_preprocessed_data("config/run-files/mayank-merged.yaml")
var_data = vec(var(data_from_dataframe(norm_counts),dims=2 ))

#histogram to see where coding and noncoding is represented in variance data
bins = Int.(floor.(var_data.*10)).+1
barplot(bins,repeat([1],length(var_data)),stack=bins,color=replace(norm_counts.transcript_type,"coding"=>1,"noncoding"=>2))




