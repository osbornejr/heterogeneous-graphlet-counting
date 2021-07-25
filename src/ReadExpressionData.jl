using DataFrames 
import CSV
    
function read_count_data(data_dir::String;method::String)
    ##method can be "expected_count","TPM", "FPKM" or "IsoPct"
    samples=CSV.read.(readdir(data_dir,join=true))
    sample_names=replace.(readdir(data_dir),"_RSEM.isoforms.results"=>"").*"_data"
    transcript_names=select(samples[1],1)
    counts=rename!(hcat(transcript_names,select.(samples,Symbol(method))...,makeunique=true),["transcript_id";sample_names])
end
    
## This function expects count data generated by the above read_count_data function, as well as file path to a list of hit 
## transcripts corresponding to transcripts in the original data. TODO this process probably should occur further up the line
## (pre-quantification even).
function filter_count_data(filter_list_path,count_data)
    filter_list=Array(DataFrame!(CSV.File(filter_list_path,header=false)))
    filtered_counts=count_data[findall(in(filter_list),count_data.transcript_id),:]
    return filtered_counts
end
