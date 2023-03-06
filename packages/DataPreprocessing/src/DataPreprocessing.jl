module DataPreprocessing
export data_from_dataframe

include("Normalisation.jl")
include("ReadExpressionData.jl")

function data_from_dataframe(df::DataFrame,identifier::String)
    #transform the numerical data columns of a DataFrame into an array. The identifier is a string that is common and uniquely contained in the names of the desired columns 
    data = Array(select(df,filter(x->occursin(identifier,x),names(df))))
end

function clean_raw_counts(raw_counts::DataFrame,expression_cutoff::Int)

    ##Deduplicate-- there may be multiple entries for the same transcript, we need to select only one of these

    ## Clean - remove transcripts with total counts across all samples less than Cut
    ##plot before cut
    #  histogram(DataFrame([log2.(vec(sum(raw_data,dims=2))),raw_counts[!,:transcript_type]],[:sum,:transcript_type]),:sum,:transcript_type,"$(outdir)/raw_data_histogram.svg",xaxis =" sum of expression (log2 adjusted)")
    raw_data = data_from_dataframe(raw_counts,"data")
    clean_counts=raw_counts[vec(sum(raw_data,dims = 2 ).>= expression_cutoff),:]

    ##plot after cut
    # histogram(DataFrame([log2.(vec(sum(clean_data,dims=2))),clean_counts[!,:transcript_type]],[:sum,:transcript_type]),:sum,:transcript_type,"$(outdir)/clean_data_cut_histogram.svg",xaxis =" sum of expression (log2 adjusted)")
    #boxplot(raw_counts,"raw_data_cleaned_boxplot.svg")
    return clean_counts
end



function normalise_clean_counts(clean_counts::DataFrame,norm_method::String)
### Normalisation

    clean_data = data_from_dataframe(clean_counts,"data")
    norm_data = library_size_normalisation(clean_data,norm_method)
    norm_counts = copy(clean_counts)
    norm_counts[:,findall(x->occursin("data",x),names(norm_counts))] = norm_data
    return norm_counts
end

function sample_norm_counts(norm_counts::DataFrame,variance_percent::Float64)
##Sampling for most variable transcripts
#add variance column to normalised data
    norm_data = data_from_dataframe(norm_counts,"data")
    variance = vec(var(norm_data, dims=2))
    norm_counts.variance = variance
    sample_counts_noncoding = sort(norm_counts[norm_counts[!,:transcript_type].=="noncoding",:],:variance)[Int(round(end*(1-variance_percent))):end,:]
    sample_counts_coding = sort(norm_counts[norm_counts[!,:transcript_type].=="coding",:],:variance)[Int(round(end*(1-variance_percent))):end,:]
    sample_counts = outerjoin(sample_counts_noncoding,sample_counts_coding,on = names(norm_counts))
    return sample_counts
end

function preprocess_raw_counts(raw_counts::DataFrame,expression_cutoff::Int,norm_method::String,variance_percent::Float64)
    #all in one function to get data processed for network construction
    return sample_norm_counts(normalise_clean_counts(clean_raw_counts(raw_counts,expression_cutoff),norm_method),variance_percent)
end
end # module
