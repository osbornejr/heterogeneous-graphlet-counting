using CSV
using DataFrames

function decide_blastn_matches(input_file::String)
    df_in = CSV.read(input_file,DataFrame,stringtype=String)
    df_out = df_in#[:,[:QueryID,:SeqID,:RNAType]]
    ##add column to denote match type
    df_out.MatchType = fill("no match",size(df_out,1))
    ##allow missing values in query column
    df_out.QueryID = convert(Vector{Union{eltype(df_out.QueryID), Missing}}, df_out.QueryID)
    relative_size = ""
    for (i,row) in enumerate(eachrow(df_in))
        if (!ismissing(row.SeqID))
            n,d = parse.(Int,split(row.TotalIdentities,"/"))
            ##if a potential match was found by blastn
            if (row.QueryLen*0.8) > row.SeqLen
                ## the query is bigger than sequence
                relative_size = "super"
                if (n > d*0.8) & (d > row.SeqLen*0.8)
                   ##there is a significant part of query matching sequence, and             at least half of those matches are not gapped).
                    df_out[i,:MatchType] = "$relative_size match" 
                else
                    ## either too many gaps, or not enough matching, so we disregard match
                    df_out[i,:MatchType] = "minor association"
                end
            elseif (row.QueryLen*1.2) < row.SeqLen
                ## the query is smaller than sequence
                relative_size = "sub"
                if (n > d*0.8) & (d > row.QueryLen*0.8)
                   ##there is a significant part of query matching sequence, and             at least half of those matches are not gapped).
                    df_out[i,:MatchType] = "$relative_size match" 
                else
                    ## either too many gaps, or not enough matching, so we disregard match
                    df_out[i,:MatchType] = "minor association"
                end
            else
                ##the query and sequence are of similar size 
                if (n > d*0.8) & (d > row.QueryLen*0.8)
                   ##there is a significant part of query matching sequence, and             at least half of those matches are not gapped).
                    df_out[i,:MatchType] = "match" 
                else
                    ## either too many gaps, or not enough matching, so we disregard match
                    df_out[i,:MatchType] = "minor association"
                end
            end
        end
        ## no match found by blastn
    end
    return df_out
end

function get_rnatype_from_description(df_out)
    for (i,row) in enumerate(eachrow(df_out))
        ## other housekeeping. check description for missing terms to put in RNAType field
        if (ismissing(row.RNAType))
            if occursin("DNA",row.Description)
                df_out[i,:RNAType] = "DNA"
            end
            if occursin("cDNA",row.Description)
                df_out[i,:RNAType] = "cDNA"
            end
            if occursin("mitochond",row.Description)
                df_out[i,:RNAType] = "mitoc"
            end
            if occursin("satellite",row.Description)
                df_out[i,:RNAType] = "satellite"
            end
            if occursin(" transposon ",row.Description)
                df_out[i,:RNAType] = "transposon"
            end
            if occursin("gypsy like",row.Description)
                df_out[i,:RNAType] = "transposon"
            end
            if occursin("retrotransposon",row.Description)
                df_out[i,:RNAType] = "retrotransposon"
            end
            if occursin("Retrotransposon",row.Description)
                df_out[i,:RNAType] = "RT gene"
            end
            if occursin("RT gene",row.Description)
                df_out[i,:RNAType] = "retrotransposon"
            end
            if occursin("pseudogene",row.Description)
                df_out[i,:RNAType] = "pseudogene"
            end
            if occursin("genomic sequence",row.Description)
                df_out[i,:RNAType] = "gene"
            end
            if occursin("gene complete cds",row.Description)
                df_out[i,:RNAType] = "gene"
            end
            if occursin("protein gene",row.Description)
                df_out[i,:RNAType] = "gene"
            end
            if occursin("complete genome",row.Description)
                df_out[i,:RNAType] = "genome"
            end
            if occursin("whole genome",row.Description)
                df_out[i,:RNAType] = "genome"
            end
            if occursin("complete sequence",row.Description)
                df_out[i,:RNAType] = "genome"
            end
            if occursin("COMPLETE SEQUENCE",row.Description)
                df_out[i,:RNAType] = "genome"
            end
            if occursin("chromosome",row.Description)
                df_out[i,:RNAType] = "chromosome"
            end
            if occursin("ribosomal RNA", row.Description)
                df_out[i,:RNAType] = "rrna"
            end
            if occursin("misc_RNA",row.Description)
                df_out[i,:RNAType] = "misc_RNA"
            end
            if occursin("long non-coding RNA",row.Description)
                df_out[i,:RNAType] = "lncrna"
            end
        end
    end
    return df_out
end

function sort_by_QueryID(df_a::DataFrame,df_b::DataFrame)
    return df_b[indexin(df_a.QueryID,df_b.QueryID),:]
end

function merge_np_and_top_matches(df_np,df_top)
    match_df = df_np
    ##get transcripts where there was a match in cicer_top but not cicer_np
    idx = findall((df_np.MatchType.=="match").-(df_top.MatchType.=="match").==-1)
    match_df[idx,:] = df_top[idx,:]
    ##now get transcripts where there was a super or sub match in cicer top but no or minor match in cicer np
    idx = findall(((df_np.MatchType.=="match").+(df_np.MatchType.=="super match").+ (df_np.MatchType.=="sub match")).-((df_top.MatchType.=="super match") .+ (df_top.MatchType.=="sub match")).==-1)
    match_df[idx,:] = df_top[idx,:]
    return match_df
end

function decide_and_merge_blastn_matches(db::String;path="data/mayank-de-novo")
    cicer_np_match_df = decide_blastn_matches("$(path)/$(db)_cicer_nonpredicted_filtered_output_final.txt")
    cicer_top_match_df = decide_blastn_matches("$(path)/$(db)_cicer_top_filtered_output_final.txt")
    any_np_match_df = decide_blastn_matches("$(path)/$(db)_any_nonpredicted_filtered_output_final.txt")
    any_top_match_df = decide_blastn_matches("$(path)/$(db)_any_top_filtered_output_final.txt")
    
    ##ensure all are sorted the same as first
    cicer_top_match_df = sort_by_QueryID(cicer_np_match_df,cicer_top_match_df)
    any_np_match_df = sort_by_QueryID(cicer_np_match_df,any_np_match_df)
    any_top_match_df = sort_by_QueryID(cicer_np_match_df,any_top_match_df)
    
    
    #start merging in all. by default we accept the cicer np at first, and improve it with the other matches if possible
    ###MERGE CICER NP <- CICER TOP
    cicer_match_df = merge_np_and_top_matches(cicer_np_match_df,cicer_top_match_df)
    ##MERGE ANY NP <- ANY TOP
    any_match_df = merge_np_and_top_matches(any_np_match_df,any_top_match_df)
    
    ##MERGE CICER <- ANY
    #can use same logic as np <- top here
    match_df = merge_np_and_top_matches(cicer_match_df,any_match_df)
    return match_df
end

##wrap load in of all original output dfs here
#begin
#    ###Load all 8 original blastn outputs for checking/comparsions
#    rs_cicer_np_df = CSV.read("data/mayank-de-novo/ref-seq_cicer_nonpredicted_filtered_output_final.txt",DataFrame,stringtype=String)
#    rs_cicer_top_df = sort_by_QueryID(rs_cicer_np_df,CSV.read("data/mayank-de-novo/ref-seq_cicer_top_filtered_output_final.txt",DataFrame,stringtype=String))
#    rs_any_np_df = sort_by_QueryID(rs_cicer_np_df,CSV.read("data/mayank-de-novo/ref-seq_any_nonpredicted_filtered_output_final.txt",DataFrame,stringtype=String))
#    rs_any_top_df = sort_by_QueryID(rs_cicer_np_df,CSV.read("data/mayank-de-novo/ref-seq_any_top_filtered_output_final.txt",DataFrame,stringtype=String))
#    nt_cicer_np_df = sort_by_QueryID(rs_cicer_np_df,CSV.read("data/mayank-de-novo/core-nt_cicer_nonpredicted_filtered_output_final.txt",DataFrame,stringtype=String))
#    nt_cicer_top_df = sort_by_QueryID(rs_cicer_np_df,CSV.read("data/mayank-de-novo/core-nt_cicer_top_filtered_output_final.txt",DataFrame,stringtype=String))
#    nt_any_np_df = sort_by_QueryID(rs_cicer_np_df,CSV.read("data/mayank-de-novo/core-nt_any_nonpredicted_filtered_output_final.txt",DataFrame,stringtype=String))
#    nt_any_top_df = sort_by_QueryID(rs_cicer_np_df,CSV.read("data/mayank-de-novo/core-nt_any_top_filtered_output_final.txt",DataFrame,stringtype=String))
#end

@info "merging blastn hits from ref-seq database..."
rs_match_df = decide_and_merge_blastn_matches("ref-seq")

@info "merging blastn hits from core-nt database..."
nt_match_df = decide_and_merge_blastn_matches("core-nt")

##sort nt on ref-seq
nt_match_df = sort_by_QueryID(rs_match_df,nt_match_df)


##MERGE REF_SEQ <- CORE_NT 
@info "merging core-nt hits into ref-seq hits..."
match_df = rs_match_df
### we will definitely take any core-nt matches and sub/super matches over a ref-seq minor or no association:
idx = findall((((rs_match_df.MatchType.=="sub match").+(rs_match_df.MatchType.=="super match").+(rs_match_df.MatchType.=="match")).-(nt_match_df.MatchType.=="match")).==-1)
match_df[idx,:] = nt_match_df[idx,:]
idx = findall((((rs_match_df.MatchType.=="sub match").+(rs_match_df.MatchType.=="super match").+(rs_match_df.MatchType.=="match")).-(nt_match_df.MatchType.=="super match")).==-1)
match_df[idx,:] = nt_match_df[idx,:]
idx = findall((((rs_match_df.MatchType.=="sub match").+(rs_match_df.MatchType.=="super match").+(rs_match_df.MatchType.=="match")).-(nt_match_df.MatchType.=="sub match")).==-1)
match_df[idx,:] = nt_match_df[idx,:]

##for replacing rs sub/super matches with nt matches, we review first to see what is replacing what
idx = findall(((rs_match_df.MatchType.=="super match").+(nt_match_df.MatchType.=="match")).==2)
##observe (for instance)
rs_any_np_df[idx,:]
nt_any_np_df[idx,:]
##by default we include the nt matches, as they should be good, and we will keep the rs and nt match outputs around as well if we do need purely "rna" matches.
match_df[idx,:] = nt_match_df[idx,:]
idx = findall(((rs_match_df.MatchType.=="sub match").+(nt_match_df.MatchType.=="match")).==2)
##observe (for instance)
rs_any_np_df[idx,:]
nt_any_np_df[idx,:]
##by default we include the nt matches, as they should be good, and we will keep the rs and nt match outputs around as well if we do need purely "rna" matches.
match_df[idx,:] = nt_match_df[idx,:]

##nt minor over rs none
##lastly, we want to take any minor associations in nt that are no matches in refseq.
idx = findall(((rs_match_df.MatchType.=="no match").+(nt_match_df.MatchType.=="minor association")).==2)
match_df[idx,:] = nt_match_df[idx,:]

@info "searching for more rnatypes in description field..."
match_df = get_rnatype_from_description(match_df)

##manually inspect the last few explicit matches that haven't got an explicit RNAtype 
filter(:MatchType=>==("match"),match_df[(ismissing.(match_df.RNAType).+ismissing.(match_df.SeqID)).==1,:])



