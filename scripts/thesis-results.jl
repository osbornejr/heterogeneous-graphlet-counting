using Pkg
Pkg.activate("results")
#dev packages
using ProjectFunctions
using DataPreprocessing
using NetworkConstruction
using GraphletCounting
using GraphletAnalysis
using Results

##external packages
using StatsBase
#using GLMakie
using CairoMakie


##set output base dir of thesis
thesis_dir = "/Users/osbornejr/git/phd/thesis/"
thesis_version_dir = "final-draft"
network_chapter_dir = "network-analysis"


#set colours
colour_map = Dict("noncoding"=>:cyan,"coding"=>:purple)

## set experiment
experiment = "mayank-merged"
experiment = "GSE68559_sub"

## Chickpea salt stress
#run = "mayank-merged-1400-network"
#run = "mayank-unmerged"
#run = "mayank-merged-altered-network"

#PCIT method
#run = "mayank-merged-small-pruned"
#run = "mayank-merged-large-pruned"

# Ridge partial correlation method
#run = "mayank-merged"

##Human smoker
#run = "GSE68559_sub"
run = "Milestone-3-network"

##or batch wise
#batch = nothing 
batch = [("GSE68559_sub","Milestone-3-network"),("mayank-merged","mayank-merged-large-pruned")]

if batch == nothing
    #construct batch just for one experiment
    batch = [(experiment,run)]
end

for er in batch
    experiment = er[1]
    run = er[2]
    ## Chickpea salt stress
    #run = "mayank-merged-1400-network"
    #run = "mayank-unmerged"
    #run = "mayank-merged-altered-network"
    
    #PCIT method
    #run = "mayank-merged-small-pruned"
    #run = "mayank-merged-large-pruned"
    
    ## Ridge partial correlation method
    #run = "mayank-merged"
    
    
    #run = "GSE68559_sub"         
    
    ##load preprocessing data
    raw_counts,round_counts,vst_counts,clean_counts,norm_counts,processed_counts = get_preprocessed_data("config/run-files/$(experiment)/$(run).yaml")
    #plot variance histogram
    Results.variance_histogram(norm_counts)
    
    
    ##load network info 
    components,adj_matrix,network_counts,vertexlist,edgelist = get_network_construction()
    #visualise network
    fig = Results.plot_network(adj_matrix,vertex_colors = replace(vertexlist,collect(colour_map)...))
    fig = Results.plot_network(adj_matrix,vertex_colors = replace(vertexlist,collect(colour_map)...),layout="Grouped",groupby = replace(vertexlist,"coding"=>0,"noncoding"=>1))
    #view typed degree distributions
    f = Results.typed_degree_distribution(vertexlist,edgelist)
    
    #for each component, dd and network density
    dds = map(x->sum(adj_matrix,dims=1)[x],components)
    
    #Results.add_to_fig(fig)
    #
    
    #finding high degree nodes
    #degree distribution
    dd = sum.(values.(GraphletCounting.typed_degree_distribution(vertexlist,edgelist)))
    #get ids for large component
    #largest component
    #l_comp = components[findall(==(maximum(length.(components))),length.(components))...]
    id_cut = 100
    #high_ids = l_comp[dd.>id_cut]
    #ids = l_comp[dd.<id_cut]
    #shouldn't need high id distinction now
    high_ids = (1:size(adj_matrix,1))[dd.>id_cut]
    ids = (1:size(adj_matrix,1))[dd.<id_cut]
    #compare processed_data for high degree to rest of proto-network
    pd = data_from_dataframe(processed_counts)
    
    ## at this stage lets rearrange for the visualisation to have the sample columns in order (baked in order is 10,11,12,1,2,3,4,5,6,7,8,9; this translates to ICCV2-Stressed,JG11-Control,JG11-Stressed,ICCV2-Control. we will move the first three columns to the back to get each strain next to its condition pair.
    #pd = pd[:,[4,5,6,7,8,9,10,11,12,1,2,3]]
    #sample_names = ["JG11-Control-R1","JG11-Control-R2","JG11-Control-R3","JG11-Stress-R1","JG11-Stress-R2","JG11-Stress-R3","ICCV2-Control-R1","ICCV2-Control-R2","ICCV2-Control-R3","ICCV2-Stress-R1","ICCV2-Stress-R2","ICCV2-Stress-R3"]
    
    f2 = Figure()
    ax1 = Axis(f2[1,1],xticks=((2:3:12),["Control","Stress","Control","Stress"]),title="Expression counts: high degree transcripts")
    ##add extra layer of xtick labels
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
    cor_threshold = 0.999
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
    binned_pd = (pd.>2)+(pd.>7)+(pd.>12)
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
    #pd = pd[:,[4,5,6,7,8,9,10,11,12,1,2,3]]
    
    #Chickpea
    #strain = ["JG11 (salt tolerant)"=>[1,2,3,4,5,6],"ICCV2 (salt sensitive)"=>[7,8,9,10,11,12]]
    #treatment = ["Control"=>[1,2,3,7,8,9],"Salt Stress"=>[4,5,6,10,11,12]]
    # Human
    if occursin("GSE68559",experiment)
      condition = "treatment" 
       elseif occursin("mayank",experiment)
           condition = "strain"
       else
           throw(ArgumentError("need to provide a condition for current experiment"))
    end
    
    
    matching_condition_pattern = DataPreprocessing.high_contrast_transcripts(pd,collect(values(params["conditions"][condition]))[1],collect(values(params["conditions"][condition]))[2],2.0,strictness="n-1")
    
    
    f4 = Figure()
    ax8 = Axis(f4[1,1],xticks=((2:3:12),["Control","Stress","Control","Stress"]),title="Expression profiles of selected transcripts")
    Axis(f4[1,1],xticks=((2:6:12),["JG11 (salt tolerant)","ICCV2 (salt sensitive)"]),xticklabelpad=25)
    ## if we want all epxression profiles shown, with the problematic onees at the top
    to_view = pd[matching_condition_pattern,:]'
    #to_view = pd[setdiff(1:size(pd,1),matching_condition_pattern),:]'
    #to_view = pd'
    pd = data_from_dataframe(processed_counts)
    to_view = pd'
    #to_view = binned_pd
    heatmap!(ax8,to_view)
    ##generalise-- finding high contrast counts. condition(s) supplied in param config?
    
    
    ## Analysis
    #
    #
    #Communities
    comm_df = ProjectFunctions.community_analysis(network_counts,adj_matrix);
    vertex_colors =string.(comm_df.color);
    c_comp = components[findall(==(maximum(length.(components))),length.(components))...]
    #for some reason communities aren't detected on largest component TODO fix/do community detection on each component
    #HACK for now, set it manually
    if length(components)>1
        c_comp = components[3]
    end
    comm_fig = Results.plot_network(adj_matrix[c_comp,c_comp],vertex_colors = vertex_colors,groupby=comm_df.group)
    
    #WGCNA
    #wgcna_network,wgcna_comm = ProjectFunctions.get_wgcna();
    #unweighted_wgcna= (wgcna_network.>0.22)
    #wgcna_components = NetworkConstruction.network_components(unweighted_wgcna)
    ##largest = findmax(length.(wgcna_components))[2]
    #wgcna_adj = unweighted_wgcna
    ##wgcna_comp_comms = wgcna_comm[wgcna_components[largest],:]
    #vertex_colors =string.(wgcna_comm.color);
    #wgcna_fig = Results.plot_network(wgcna_adj,vertex_colors = vertex_colors)

    

    ##Network vis for thesis
    save("$(thesis_dir)/$(thesis_version_dir)/$(network_chapter_dir)/figs/$(experiment)_network.png",comm_fig,px_per_unit =4.0)




    #
    #Typed representations
    graphlet_counts,timer = get_graphlet_counts()
    t_r_output = typed_representations(graphlet_counts,timer,vertexlist,edgelist)
    t_r_result = Results.typed_representation_results(t_r_output,colour_mapping = colour_map)
    save("$(thesis_dir)/$(thesis_version_dir)/$(network_chapter_dir)/figs/$(experiment)_t_r_results.pdf",t_r_result)
end 
