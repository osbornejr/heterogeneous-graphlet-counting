module DataPreprocessing
export data_from_dataframe

include("Normalisation.jl")
include("ReadExpressionData.jl")
include("Visualisation.jl")

function data_from_dataframe(df::DataFrame,identifier::String="data")
    #transform the numerical data columns of a DataFrame into an array. The identifier is a string that is common and uniquely contained in the names of the desired columns 
    data = Array(select(df,filter(x->occursin(identifier,x),names(df))))
end
export data_from_dataframe

function round_raw_counts(raw_counts::DataFrame,sig::Int)
    raw_data = data_from_dataframe(raw_counts,"data")
    round_data = round.(raw_data,digits=5)
    round_counts = copy(raw_counts)
    round_counts[:,findall(x->occursin("data",x),names(round_counts))] = round_data
    return round_counts
end

function log_counts(counts::DataFrame,pseudocount::Float64=1.0)
    data = data_from_dataframe(counts,"data")
    log_data = log2.(data.+pseudocount)
    log_counts = copy(counts)
    log_counts[:,findall(x->occursin("data",x),names(log_counts))] = log_data
    return log_counts
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
    elseif (method == "old_method")
        ##adding here for testing purposes. it uses the old method which only requires one input, expression_cutoff which should instead be an integer-- the cut_percent parameter is still labelled expression cutoff in configs due to this-- for testing we keep that match and ignore minreq. All transcripts with total counts below this figure will be removed. It was usually applied to raw, not rounded counts, but that shouldn't make too much difference.
        clean_counts = clean_raw_counts(round_counts,Integer(cut_percent))
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

function sample_norm_counts(norm_counts::DataFrame,variance_cutoff::Float64;method::String="cutoff",maintain_ratio=true)

##Sampling for most variable transcripts
#add variance column to normalised data
    norm_data = data_from_dataframe(norm_counts,"data")
    variance = vec(var(norm_data, dims=2))
    if (method == "percent")
        if (variance_cutoff>1) | (variance_cutoff<0)
            throw(ArgumentError("'percent' method requires a variance percentage input"))
        end
        variance_percent = variance_cutoff
        norm_counts.variance = variance
        #TODO generalise to any types/number of types
        if (maintain_ratio == true)
            sample_counts_noncoding = sort(norm_counts[norm_counts[!,:transcript_type].=="noncoding",:],:variance)[Int(round(end*(1-variance_percent))):end,:]
            sample_counts_coding = sort(norm_counts[norm_counts[!,:transcript_type].=="coding",:],:variance)[Int(round(end*(1-variance_percent))):end,:]
            sample_counts = outerjoin(sample_counts_noncoding,sample_counts_coding,on = names(norm_counts))
        else
            sample_counts = sort(norm_counts,:variance)[Int(round(end*(1-variance_percent))):end,:]
        end
    elseif (method == "cutoff")
        # new, simpler method just setting variance cutoff
        sample_counts = norm_counts[variance.>variance_cutoff,:]
    else
        throw(ArgumentError("method must be either 'cutoff' or 'percent'"))
    end

    return sample_counts
end

end # module
