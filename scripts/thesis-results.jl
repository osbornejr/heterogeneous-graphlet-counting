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
fig_dir = "$(thesis_dir)/$(thesis_version_dir)/$(network_chapter_dir)/figs/"

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
    #data set for step before variance cut (after normalisation) 
    pd = data_from_dataframe(norm_counts)
    
    #set condition based on experiment
    if occursin("GSE68559",experiment)
      condition = "treatment" 
       elseif occursin("mayank",experiment)
           condition = "strain"
       else
           throw(ArgumentError("need to provide a condition for current experiment"))
    end
    
    
    ##define rev dictionary here for treatment values (used for figure)
    rev_dict = Dict()
    for (k,v) in params["conditions"]["treatment"]
        for x in v
            push!(rev_dict,x=>k)
        end
    end

    # we also need to sort samples based on their treatment type, within each of the chosen condition.
    display_sample_order = []
    display_condition_order = []
    display_condition_ticks = []
    display_treatment_order = []
    display_treatment_ticks = []
    for (k,v) in params["conditions"][condition]
        treatment_match = [(x,rev_dict[x]) for x in v]
        sort!(treatment_match, by = x->x[2])
        push!(display_sample_order,first.(treatment_match)...)
        push!(display_condition_order,k)
        push!(display_condition_ticks,length(v)/2+.5)
        for (k,v)in countmap(last.(treatment_match))
            push!(display_treatment_ticks,v/2+.5)
            push!(display_treatment_order,k)
        end
    end


    ##if chosen condition has multiple features, we need to run for each condition against all others (i.e. region 1 vs region 2-10, for regions 1-10)
    ## two set condition is actually exception (we don't need to compare each set against the other twice)
    if length(params["conditions"][condition]) == 2
        condition_a =collect(values(params["conditions"][condition]))[1]
        condition_b =collect(values(params["conditions"][condition]))[2]
   ##find which transcripts follow the chosen pattern 
    matching_condition_pattern = DataPreprocessing.high_contrast_transcripts(pd,condition_a,condition_b,2.0,strictness="n-1")
    fig = Figure()
    # we apply treatment level axis first
    ax = Axis(fig[1,1],xticks=(display_treatment_ticks,display_treatment_order),title="Expression profiles of selected transcripts")
    Axis(f4[1,1],xticks=(display_condition_ticks,display_condition_order),xticklabelpad=25)
    heatmap!(ax,pd[:,display_sample_order])
    else
    
        for (cond_name,cond_samples) in params["conditions"][condition]
            condition_a = cond_samples
            ## other condition is all other samples
            condition_b = setdiff(vcat(values(params["conditions"][condition])...),condition_a)
            ##find which transcripts follow the chosen pattern 
            matching_condition_pattern = DataPreprocessing.high_contrast_transcripts(pd,condition_a,condition_b,2.0,strictness="n-1")
            fig = Figure()
            # we apply treatment level axis first
            ax = Axis(fig[1,1],xticks=(display_treatment_ticks,display_treatment_order),title="Expression profiles of selected transcripts")
            Axis(f4[1,1],xticks=(display_condition_ticks,display_condition_order),xticklabelpad=25)
            heatmap!(ax,pd[:,display_sample_order])
        end
    end

    #TODO
## - adjust to allow for multifeature conditions (as in human data)
#  - cycle through all conditions (should also happen in run)
## - set labels on fig based on actual feature info
## - set for publication quality
#       - as pdf
#       - right background
#       - latex font

    
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
    save("$(fig_dir)/$(experiment)_network.png",comm_fig,px_per_unit =4.0)



    #
    #Typed representations
    graphlet_counts,timer = get_graphlet_counts()
    t_r_output,merged_summaries = typed_representations(graphlet_counts,timer,vertexlist,edgelist)
    ## table of over/under represented
    t_r_result = Results.typed_representation_results(t_r_output,colour_mapping = colour_map)
    save("$(fig_dir)/$(experiment)_t_r_results.pdf",t_r_result)
    #merged boxplots
    NetworkConstruction.tex_merged_boxplot(merged_summaries,"$(fig_dir)/$(experiment)_merged_boxplot.tex","input",ylabel = "log value")
end 
