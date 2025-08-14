using Pkg
Pkg.activate(".")
using ProjectFunctions
using GraphletCounting
using StatsBase
using Results
using GLMakie
##set experiment
experiment = "mayank-merged-1400-network"
#experiment = "mayank-unmerged"
#experiment = "mayank-merged-altered-network"
#experiment = "mayank-merged"
#experiment = "GSE68559_sub"

##load preprocessing data
raw_counts,round_counts,vst_counts,clean_counts,norm_counts,processed_counts = get_preprocessed_data("config/run-files/$(experiment).yaml")
#plot variance histogram
Results.variance_histogram(norm_counts)


##load network info 
components,adj_matrix,network_counts,vertexlist,edgelist = get_network_construction()
#visualise network
fig = Results.plot_network(adj_matrix,vertex_colors = replace(vertexlist,"noncoding"=>:blue,"coding"=>:purple))
#view typed degree distributions
f = Results.typed_degree_distribution(vertexlist,edgelist)

#Results.add_to_fig(fig)
#


#finding high degree nodes
#degree distribution
dd = sum.(values.(GraphletCounting.typed_degree_distribution(vertexlist,edgelist)))
#get ids for large component
#largest component
l_comp = components[findall(==(maximum(length.(components))),length.(components))...]
id_cut = 100
high_ids = l_comp[dd.>id_cut]
ids = l_comp[dd.<id_cut]
#compare processed_data for high degree to rest of proto-network
pd = data_from_dataframe(processed_counts)

## at this stage lets rearrange for the visualisation to have the sample columns in order (baked in order is 10,11,12,1,2,3,4,5,6,7,8,9; this translates to ICCV2-Stressed,JG11-Control,JG11-Stressed,ICCV2-Control. we will move the first three columns to the back to get each strain next to its condition pair.
pd = pd[:,[4,5,6,7,8,9,10,11,12,1,2,3]]
sample_names = ["JG11-Control-R1","JG11-Control-R2","JG11-Control-R3","JG11-Stress-R1","JG11-Stress-R2","JG11-Stress-R3","ICCV2-Control-R1","ICCV2-Control-R2","ICCV2-Control-R3","ICCV2-Stress-R1","ICCV2-Stress-R2","ICCV2-Stress-R3"]

f2 = Figure()
ax1 = Axis(f2[1,1],xticks=((2:3:12),["Control","Stress","Control","Stress"]),title="Expression counts: high degree transcripts")
#add extra layer of xtick labels
Axis(f2[1,1],xticks=((2:6:12),["JG11 (salt tolerant)","ICCV2 (salt sensitive)"]),xticklabelpad=25)
heatmap!(ax1,pd[high_ids,:]')
ax2 = Axis(f2[1,2],xticks=((2:3:12),["Control","Stress","Control","Stress"]),title="Expression counts: other transcripts")
Axis(f2[1,2],xticks=((2:6:12),["JG11 (salt tolerant)","ICCV2 (salt sensitive)"]),xticklabelpad=25)
heatmap!(ax2,pd[ids,:]')
f2

# it is clear that PCIT is not working as it should, and it appears the issue is the parital correlation step. If we try to get partialcor matrix manually in julia if also fails, as the underlying cor matrix is not PSD and thus not invertible. The underlying issue seems to be too many highly correlated transcripts. Lets identify them here.

#take correlation matrix
cor_pd = cor(pd')
#find all transcript 9pairs with cor greater than threshold
cor_threshold = 0.99
matching_pairs = filter(t->t[1]!=t[2],Tuple.(findall(>(cor_threshold),cor_pd)))
high_cors = unique(vcat(first.(matching_pairs),last.(matching_pairs)))
ax3 = Axis(f2[2,1],xticks=((2:3:12),["Control","Stress","Control","Stress"]),title="Expression counts: highly correlated (>$(cor_threshold)) transcripts")
Axis(f2[2,1],xticks=((2:6:12),["JG11 (salt tolerant)","ICCV2 (salt sensitive)"]),xticklabelpad=25)
heatmap!(ax3,pd[high_cors,:]')

##observe those transcripts that have high cor interactions but not high degree 
non_high_deg_high_cors = high_cors[.!in.(high_cors,Ref(intersect(high_ids,high_cors)))]
ax4 = Axis(f2[2,2],xticks=((2:3:12),["Control","Stress","Control","Stress"]),title="Expression counts: highly correlated but not high degree transcripts")
Axis(f2[2,2],xticks=((2:6:12),["JG11 (salt tolerant)","ICCV2 (salt sensitive)"]),xticklabelpad=25)
heatmap!(ax4,pd[non_high_deg_high_cors,:]')

non_high_cor_high_deg = high_ids[map(x->!in(x,intersect(high_ids,high_cors)),high_ids)]
ax5 = Axis(f2[3,1],xticks=((2:3:12),["Control","Stress","Control","Stress"]),title="Expression counts: high degree but not high cor transcripts")
Axis(f2[3,1],xticks=((2:6:12),["JG11 (salt tolerant)","ICCV2 (salt sensitive)"]),xticklabelpad=25)
heatmap!(ax5,pd[non_high_cor_high_deg,:]')

##find those transcripts that have an expression pattern that is on or off depending on strain (with a tolerance of 1 sample not following pattern)
#matching_strain_pattern = collect(1:size(pd)[1])[abs.(sum((pd.==0.0)[:,1:6],dims=2)-sum((pd.==0.0)[:,7:12],dims=2)).>4]
matching_JG_strain_pattern = collect(1:size(pd)[1])[sum((pd.==0.0)[:,1:6],dims=2)-sum((pd.==0.0)[:,7:12],dims=2).>4]
matching_ICCV_strain_pattern = collect(1:size(pd)[1])[sum((pd.==0.0)[:,1:6],dims=2)-sum((pd.==0.0)[:,7:12],dims=2).<-4]
matching_strain_pattern = vcat(matching_JG_strain_pattern,matching_ICCV_strain_pattern)
ax6 = Axis(f2[3,2],xticks=((2:3:12),["Control","Stress","Control","Stress"]),title="Expression counts: expression profiles matching strain samples (with 1 sample tolerance)")
Axis(f2[3,2],xticks=((2:6:12),["JG11 (salt tolerant)","ICCV2 (salt sensitive)"]),xticklabelpad=25)
heatmap!(ax6,pd[matching_strain_pattern,:]')


### finding potentially a more agnostic way to determine which regularly occuring patterns to remove from the data
# divide the processed expression counts up into 4 bin2: x<2, 2<x<7,6<x<12, 12<x 
binned_pd = (pd.>2)#+(pd.>7)+(pd.>12)
#turn each profile into its own object
binned_profiles = map(Tuple,eachrow(binned_pd))
#sort by least to most common binned expression count profile.
sorted = sort(collect(pairs(countmap(binned_profiles))), by = x -> x[2])
#now want the transcript index so we can sort and view as heatmap
#create a Dict: profile => indices where it appears
profile_to_indices = Dict{Tuple, Vector{Int}}()
for (i, profile) in enumerate(binned_profiles)
    push!(get!(profile_to_indices, profile, Int[]), i)
end
#extract indices in sorted order
sorted_indices = [profile_to_indices[profile] for (profile, _) in sorted]
#plot processed data sorted this way as heatmap
f3 = Figure()
ax7 = Axis(f3[1,1],xticks=((2:3:12),["Control","Stress","Control","Stress"]),title="Expression counts: expression profiles binned and sorted by prevalence")
Axis(f3[1,1],xticks=((2:6:12),["JG11 (salt tolerant)","ICCV2 (salt sensitive)"]),xticklabelpad=25)
heatmap!(ax7,pd[vcat(sorted_indices...),:]')

## But actually is probably best just to select via profiles that match strains. 
#That way there is less risk of removing transcripts that are biologically interesting. 
#Let's rather set selection based on an off/on for each strain i.e. 5 or more samples on (expression>2) is on for that strain. 
#Then if we remove all transcripts that are off/on or on/off. 
#transcripts that are JG on, ICCV off
##reset back to all normalised data
pd = data_from_dataframe(norm_counts)
#rearrange for JG and ICCV order
pd = pd[:,[4,5,6,7,8,9,10,11,12,1,2,3]]

matching_JG_strain_pattern = vcat(collect(1:size(pd)[1])[(sum(pd[:,1:6].>2,dims = 2).==6).*(sum(pd[:,7:12].<2,dims=2).==6)],collect(1:size(pd)[1])[(sum(pd[:,1:6].>2,dims = 2).==5).*(sum(pd[:,7:12].<2,dims=2).>4)],collect(1:size(pd)[1])[(sum(pd[:,1:6].>2,dims = 2).>4).*(sum(pd[:,7:12].<2,dims=2).==5)])
#transcripts that are JG off, ICCV on
matching_ICCV_strain_pattern = vcat(collect(1:size(pd)[1])[(sum(pd[:,1:6].<2,dims = 2).==6).*(sum(pd[:,7:12].>2,dims=2).==6)],collect(1:size(pd)[1])[(sum(pd[:,1:6].<2,dims = 2).==5).*(sum(pd[:,7:12].>2,dims=2).>4)],collect(1:size(pd)[1])[(sum(pd[:,1:6].<2,dims = 2).>4).*(sum(pd[:,7:12].>2,dims=2).==5)])
matching_strain_pattern = vcat(matching_JG_strain_pattern,matching_ICCV_strain_pattern)

##generalise-- finding high contrast counts. condition(s) supplied in param config?
#Chickpea
strain = ["JG11 (salt tolerant)"=>[1,2,3,4,5,6],"ICCV2 (salt sensitive)"=>[7,8,9,10,11,12]]
treatment = ["Control"=>[1,2,3,7,8,9],"Salt Stress"=>[4,5,6,10,11,12]]
# Human

"""
    remove_high_contrast_transcripts(count_data,a,b,expression_cutoff;strictness,return_matches)

    Function to identify transcripts with expression profiles that follow conditions a and b closely.
    Arguments:
    - count_data: a matrix containing normalised count data, with each row representing a transcript and each column a sample.
    - a: columns corresponding to condition a.
    - b: columns corresponding to condition b.
    - expression_cutoff: the threshold for determining whether an individual count is off or on. Defaults to 0, but will be determined by the distribution of the normalised data.
    - strictness: the level at which we determine if a transcript is 'off' or 'on' for a given condition. We are identifying transcripts that are either off in condition a and on in condition b, or off in condition b and on in condition a. Defaults to 'n', i.e. that all samples in a condition must be off/on for a transcript to be off/on. The other option, 'n-1', allows a slight tolerance of one sample not bmatching the others.
    - return_matches: Whether the function returns a list of the identified transcripts alongside the pruned count data. Defaults to false.

"""
function remove_high_contrast_transcripts(count_data::Matrix{Float64},a::Vector{Int},b::Vector{Int},expression_cutoff::Float64=0;strictness::String="n",return_matches=false)
    ##set how strictly we determine an 'off' or 'on' condition-- currently must match/not match on either or all samples ("n") or all samples bar one "n-1".
    if strictness == "n"
        tol_a = length(a)
        tol_b = length(b)
    elseif strictness == "n-1"
        tol_a = length(a)-1
        tol_b = length(b)-1
    else 
        throw(ArgumentError("strictness must be either 'n' or 'n-1'."))
    end
    # identify all transcripts that are on in condition a and off in condition b 
    if tol_a == length(a) #n
        matching_a = collect(1:size(count_data)[1])[(sum(count_data[:,a].>=expression_cutoff,dims = 2).==length(a)).*(sum(count_data[:,b].<=expression_cutoff,dims=2).==length(b))]
    else #n-1
        ##maintain order, so that fully contrasted (strictness "n") transcripts are listed first, followed by either those that have exactly n-1 'on' matches in a (and any 'off' matches in b), and then by any 'on' matches in a that correspond to exactly n-1 'off' matches in b
        matching_a = vcat(collect(1:size(count_data)[1])[(sum(count_data[:,a].>=expression_cutoff,dims = 2).==length(a)).*(sum(count_data[:,b].<=expression_cutoff,dims=2).==length(b))],collect(1:size(count_data)[1])[(sum(count_data[:,a].>=expression_cutoff,dims = 2).==tol_a).*(sum(count_data[:,b].<=expression_cutoff,dims=2).>tol_b-1)],collect(1:size(count_data)[1])[(sum(count_data[:,a].>=expression_cutoff,dims = 2).>tol_a-1).*(sum(count_data[:,b].<=expression_cutoff,dims=2).==tol_b)])
    end
    # identify all transcripts that are on in condition b and off in condition a 
    if tol_b == length(b) #n
        matching_b = collect(1:size(count_data)[1])[(sum(count_data[:,b].>=expression_cutoff,dims = 2).==length(b)).*(sum(count_data[:,a].<=expression_cutoff,dims=2).==length(a))]
    else #n-1
        ##maintain order, so that fully contrasted (strictness "n") transcripts are listed first, followed by either those that have exactly n-1 'on' matches in b (and any 'off' matches in a), and then by any 'on' matches in b that correspond to exactly n-1 'off' matches in a
        matching_b = vcat(collect(1:size(count_data)[1])[(sum(count_data[:,b].>=expression_cutoff,dims = 2).==length(b)).*(sum(count_data[:,a].<=expression_cutoff,dims=2).==length(a))],collect(1:size(count_data)[1])[(sum(count_data[:,b].>=expression_cutoff,dims = 2).==tol_b).*(sum(count_data[:,a].<=expression_cutoff,dims=2).>tol_a-1)],collect(1:size(count_data)[1])[(sum(count_data[:,b].>=expression_cutoff,dims = 2).>tol_b-1).*(sum(count_data[:,a].<=expression_cutoff,dims=2).==tol_a)])
    end
    #get all transcripts that match either of above
    matching = vcat(matching_a,matching_b)
    if return_matches == true
        return [count_data[matching,:],matching]
    else
        return count_data[matching,:]
    end
end


function prune_norm_counts(norm_counts)

end
