module DataPreprocessing
export data_from_dataframe

include("Normalisation.jl")
include("ReadExpressionData.jl")

function data_from_dataframe(df::DataFrame,identifier::String)
    #transform the numerical data columns of a DataFrame into an array. The identifier is a string that is common and uniquely contained in the names of the desired columns 
    data = Array(select(df,filter(x->occursin(identifier,x),names(df))))
end

function round_raw_counts(raw_counts::DataFrame,sig::Int)
    raw_data = data_from_dataframe(raw_counts,"data")
    round_data = round.(raw_data,digits=5)
    round_counts = copy(raw_counts)
    round_counts[:,findall(x->occursin("data",x),names(round_counts))] = round_data
    return round_counts
end

function clean_round_counts(round_counts::DataFrame,cut_percent::Float64,minreq::Float64;method::String="global",output_cut::Bool=false)
    
    if (method == "global")
        ##Global cut method
        round_data = data_from_dataframe(round_counts,"data")
        n,m = size(round_data)
        ##need to get a global threshold to measure across a feature
        round_vec = sort(vec(round_data))
        ##remove zeros to determine cut
        nonzero_vec = filter(>(0),round_vec)
        # get bound on top cut_percent of ALL nonzero raw values
        cut = nonzero_vec[end-Int(round(length(nonzero_vec)*cut_percent))]
        ## get transcripts that have at least minreq of values above cut
        clean_counts = round_counts[vec(sum(round_data.>cut,dims=2).>m*minreq),:]

    elseif (method == "per-sample")
        ##per sample cut method
        round_data = data_from_dataframe(round_counts,"data")
        n,m = size(round_data)
        ##need to get threshold for each sample to check against feature
        #sort each column (sample)
        round_sorted = sort(round_data,dims=1)
        #find number of zero entries in each sample
        zs = sum(round_sorted.==0,dims=1)
        #get nonzeros for each sample
        nzs = n.-zs
        ## get row cutoff for each sample 
        cuts = n.-Int.(round.(nzs.*cut_percent))
        # get cut values for each sample
        cut = [round_sorted[cuts[i],i] for i in 1:m]
        ## get transcripts that have at least minreq of values above cut
        clean_counts = round_counts[vec(sum(round_data.>hcat([cut for i in 1:n]...)',dims=2).>m*minreq),:]
    else
        throw(ArgumentError("method must be either 'global' or 'per-sample'"))
    end
    #option to output cut or cut vector that was used to clean
    if (output_cut == true)
        return [clean_counts,cut]
    else
        return clean_counts
    end
end


function clean_raw_counts(raw_counts::DataFrame,expression_cutoff::Int)
    @info "This method is depreceated; cleaning should now happen on rounded counts using clean_round_counts"
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
