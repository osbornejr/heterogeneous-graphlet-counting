### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ 74c92bd0-ef4f-41a0-bc4b-5063ffcd12df
### set up notebook environment
begin
	import Pkg
	dev = Pkg.develop
	dev(path="packages/DataPreprocessing")
	dev(path="packages/NetworkConstruction")
	dev(path="packages/GraphletCounting")
	dev(path="packages/GraphletAnalysis")
	dev(path="packages/ProjectFunctions")
	using Revise	
	using JLD
	using JLD2
	using StatsBase
	using Gadfly,CategoricalArrays,Colors, Compose
	using DataFrames
	using DataPreprocessing
	using NetworkConstruction
	using GraphletCounting
	using GraphletAnalysis
	using ProjectFunctions
	
	cwd = ENV["JULIA_PROJECT"];
end;

# ╔═╡ 940181d4-a9b0-47e4-a13d-db2eb175e22e
### Setup input data (from human smoker GSE68559 dataset) loaded from existing cache
begin
	params = ProjectFunctions.RunParameters("GSE68559","menu1","$cwd/website",25,"upper_quartile",0.025,"pidc",0.95,"empirical_dist_zero",1000,true,false,true,true);
		

	biomart_raw_counts_file = "$cwd/output/cache/$(params.test_name)_raw_counts_biomart.jld2";
	
	raw_counts = ProjectFunctions.cache_load(biomart_raw_counts_file,"raw counts");


		
end;

# ╔═╡ 3f3e6d45-d576-4385-bea8-e55a37d34512
begin
	using CSV
	coincidents_file = "$cwd/output/cache/$(params.test_name)/cutoff/$(params.expression_cutoff)/normalisation/$(params.norm_method)/sampling/$(params.variance_percent)/similarity/$(params.coexpression)/threshold/$(params.threshold)/threshold_method/$(params.threshold_method)/analysis/graphlets/validation/coincidents/coincidents.csv";
	                        Coincidents = CSV.read(coincidents_file,DataFrame)
                        #because CSV converts the array columns to strings, we have to convert back (cost of using the easy/dirty CSV option!)
                        fix(g) = split(replace(replace(replace(replace(g,("["=>"")),("]"=>"")),("\""=>"")),(" "=>"")),",")
                        fix_int(g) = map(x->parse(Int,x),split(replace(replace(g,("["=>"")),("]"=>"")),","))
                        fix_bool(g) = BitArray(map(x->parse(Int,x),split(replace(replace(replace(g,("["=>"")),("]"=>"")),("Bool"=>"")),",")))
                        Coincidents.Vertices = fix_int.(Coincidents.Vertices)
                        Coincidents.Entrez = fix_int.(Coincidents.Entrez)
                        Coincidents.Ensembl = fix.(Coincidents.Ensembl)
                        Coincidents.Transcript_type = fix.(Coincidents.Transcript_type)
                        Coincidents.Inclusion = fix_bool.(Coincidents.Inclusion)
end;

# ╔═╡ 72b2183e-99a9-4c9b-b74b-7e7966eb4bb8
md"""
# Orbit-based biological significance
"""

# ╔═╡ 28836056-6604-4918-9f74-39bf81ad0559
md"""
## Introduction
Now that we have our heterogeneous graphlet counts and have looked at how different types are represented over the network, we want to identify the higher order interactions of particular biological relevance that might indicate the roles of non-coding transcripts with no previous functional annotation.
To do this, we will use KEGG pathway information to deterimine which biological processes are represented by the interactions in our network and infer function based on connections with known pathway transcripts.
The identified pathways will thus provide both a validation on network level that relevant biological processes have been captured by the data, and show that graphlet counting of heterogeneous networks provides an alternative method tof functional annotation that focuses much more on the local interactions between transcripts rather than the global context of gene module (community) detection which partitions the network into functional associations.        


## Candidate pathway identification 
As a first pass attempt to capture the pathway information in the network, the Entrez ids corresponding to each transcript were matched to associated KEGG pathways and the statistically significant (p<0.05) pathways over the whole network node set were selected.
This was conducted in  `R` using the `topKEGG()` function provided by the package `limma`, and the list of candiddate pathways can be seen below.
"""

# ╔═╡ 3c79ec7e-7e64-4f24-aa65-35c79167301b
md"""
## Coincident graphlets
We can use the graphlet level information in the network to deterimine these interactions beyond direct edge associations.
However, this requires a different approach to the global graphlet counts we have been working with so far.
Instead, we require for each graphlet a record of exactly which nodes were involved in each recorded graphlet relationship.
The per-edge process of the Rossi algorithm allows us to access this (at increased computational cost) within the original graphlet counting framework.
For each graphlet, we can then decide if there is association with one of our candidate pathways via the prescence of at least two known pathway nodes (we focus on __four node graphlets__ here; the level of known node representation required may change for different orders).   
A sample list of these graphlets that are **coincident** to candidate pathways are shown here:

"""

# ╔═╡ 752f5e2a-6a63-442f-9e28-90db65b6eeb1
md"""
In the above table, *Coincident_type* refers to the number of nodes in the graphlet that are known to be in the *Pathway* in question.
The other columns indidate the ids of the nodes in the graphlet in each context, save for the last column, *Inclusion*, that identifies __which__ of the nodes in the graphlet are the known pathway transcripts.
"""

# ╔═╡ f759d399-f76a-4cb5-b0d6-bf64aab1380a
sample_counts = DataPreprocessing.preprocess_raw_counts(raw_counts,params.expression_cutoff,params.norm_method,params.variance_percent);

# ╔═╡ a380f1ed-32bb-46f1-bd2c-8c633afa70d0
begin
	sim_file = "$cwd/output/cache/$(params.test_name)/cutoff/$(params.expression_cutoff)/normalisation/$(params.norm_method)/sampling/$(params.variance_percent)/similarity/$(params.coexpression)/similarity_matrix.jld2";
    similarity_matrix = cache_load(sim_file,"similarity_matrix");
end;

# ╔═╡ 6a2d3d61-cc8e-48ba-8839-705979033484
md"""
## Orbit significance
We can then look at these coincident graphlets from the perspective of each node $v$.
For each coincident graphlet that $v$ is included in, the orbital position of $v$ is recorded in relation to the associated pathway.
Each node will then have an orbit profile for each pathway, detailing how many times it is involved in association with known pathway nodes and the nature of these associations.

Importantly, we disregard the status of $v$ in regards to each pathway when considering a conincident graphlet; there must be two __other__ known pathway nodes in the graphlet.
This allows us to observe for each pathway how connected known pathway nodes are in the network.
Here we have two examples of pathways:

"""

# ╔═╡ 2abc718c-e711-4b7f-880c-38922e225dd2
md"""
Each node  that is known to be in the pathway is represented as a point for each orbit.
The height of the point indicates the number of times that node occurred in that orbit.
In the first pathway, *Glyoxylate and dicarboxylate metabolism*, there is little  signal across any of the orbits. 
In contrast, the *Alzheimer disease* pathway clearly has a lot of coincidence at the graphlet level.
We can deduce that pathways with a high graphlet connectivity in our network are those that are most active as the known pathway nodes are interacting with each other.
This means we can then identify other nodes that are also have a interaction with these known pathway nodes and infer that these nodes also have some importance in the KEGG pathway. 
Thus we narrow our scope to only those pathways that have a high connectivity across orbit categories.
The pathways we exclude may have a significant number of transcripts present in the network, but (at least at the four node level) they do not appear to be active across the sample space. 
"""


# ╔═╡ 4ee8ef83-ee2a-458c-b7ce-bf16fd4c5baa
md"""
## Removing low signal pathways
A broadbrush method to exclude these less active pathways is to identify those which have a higher percentage of **zero** orbit counts for their known pathway nodes than is observed over the whole node set.
If this occurs for a majority of four node orbits, we remove that pathway from our list.
After conducting this process, the remaining pathways of interest are

"""

# ╔═╡ 822164ee-0193-4c16-a245-b55ea8082530
md"""
## Empirical distribution of orbit counts for known pathway nodes

For these significant pathways, we want to look at the distribution of coincident orbit counts for known pathway nodes.
We can use these distributions to then compare against the coincident orbit counts of unknown nodes.
Here are some example pathway (cumulative) empirical distributions for each orbit: 


"""

# ╔═╡ 07904aca-d5da-444a-8fa5-55f6af49eb23
md"""
Note in the ecdf plot for *Fluid shear stress and atherosclerosis* above, the cumulative distribution for the `4-cycle_central` orbit is at probability 1 across the entire domain.
This indicates that there were no known pathway `4-cycle_central` counts for this pathway.
In principle, these cases are what we seek to eliminate by removing the low signal pathways above, but even after that process there may still be a few specific orbits of pathways that do not register any known pathway node coincidences.
To avoid these outliers from skewing the analysis via indicating a trivial significance to unknown nodes, we apply a filter on any such orbit-pathway pairs.
As we can see from the table below, the example shown above is the only such case in this dataset.

"""

# ╔═╡ 7f35a6f0-71e3-4448-8dd1-fe359ac73f94
md"""

## Empirical distribution comparison
We can now check every node against these empirical distributions to determine where it lies (for each orbit count) in comparison to the known pathway nodes using a threshold value $p =0.05$.
Take a node $v$ and a pathway $S$. For each orbit $H\in \mathcal{H}$, we note if the orbit count for $v$, $H_v$, is in the top $p$ of counts for $H$ across known pathway nodes of $S$.
This is done by applying the associated ecdf funtion for $S$ and $H$.
If so, we define $v$ as being $H$ significant for the pathway $S$.

If $v$ is significant for a majority of all the 11 4-node orbits in $\mathcal{H}$,
we identify $v$ as a putative member of the pathway $S$.

After applying the above process to all nodes in our network, we get the following associations:

"""

# ╔═╡ 877da2cf-fe40-4bfd-b44c-478427f3c81f
begin
		adj_file = "$cwd/output/cache/$(params.test_name)/cutoff/$(params.expression_cutoff)/normalisation/$(params.norm_method)/sampling/$(params.variance_percent)/similarity/$(params.coexpression)/threshold/$(params.threshold)/threshold_method/$(params.threshold_method)/adjacency_matrix.jld2"
                pre_adj_matrix = cache_load(adj_file,"pre-adj_matrix")
                adj_matrix = cache_load(adj_file,"adjacency_matrix")

end;

# ╔═╡ 99f1209a-3dea-4c60-8b5b-fe7da9219a1d
begin
        #Trim nodes with degree zero
        network_counts = sample_counts[vec(sum(pre_adj_matrix,dims=2).!=0),:]
        
        #maintain list of vertices in graph
        vertex_names = network_counts[!,:transcript_id]
        vertexlist = copy(network_counts[!,:transcript_type])     
        edgelist = NetworkConstruction.edgelist_from_adj(adj_matrix)
	
end;

# ╔═╡ d2f8505a-6d00-4a75-85d9-aba9bb251302


# ╔═╡ f34235b1-0e47-4f28-94db-ce3cfa598a91
begin
	
		
	    sub_Coincidents = filter(:Hom_graphlet=>x->occursin("4-",x),Coincidents);
	
end;

# ╔═╡ 08edac1b-f869-4464-828d-bcb60c64a98e
sub_Coincidents

# ╔═╡ 9bb259ec-8528-4049-80bf-5fa0d543e47c
begin
	kegg_file = "/home/osbornejr/app/output/cache/GSE68559/cutoff/25/normalisation/upper_quartile/sampling/0.025/similarity/pidc/threshold/0.95/threshold_method/empirical_dist_zero/analysis/graphlets/validation/coincidents/kegg_info.jld2"

	                    entrez_id_vector = cache_load(kegg_file,"entrez_id_vector")
                    candidates = cache_load(kegg_file,"candidates")
                    top_terms = cache_load(kegg_file,"top_terms")
                candidate_pathways = collect(keys(candidates))

end;

# ╔═╡ 8267a9ad-9551-4842-9194-5e250e79305e
top_terms

# ╔═╡ e01833de-63e0-493a-8d6c-0240c1fe1633
begin
	ermett = length(findall(x->x>0,(length.(values(GraphletAnalysis.pathways_per_node_dict([x for x in 1:length(vertexlist)],candidates))))))
	noncode_count =  countmap(vertexlist[entrez_id_vector.==0])["noncoding"]
	non_entrez = sum(entrez_id_vector.==0)
end;

# ╔═╡ 2733aa29-9e63-4b60-a7ee-0e9abcd62974
begin
md"""
Note that whilst KEGG annotations require the use of Entrez ids, there is some tension in using these gene level identifiers for transcripts.
	Several transcripts may be linked to the same gene id, and the pathways associated with that gene will thus be shared with all child transcripts, even if their individual functions differ.
$(ermett) transcripts have a match with at least one of our candidate pathways, whilst 
$(non_entrez) transcripts have no matching Entrez id; of these, the large majority
( $(noncode_count) ) are non-coding transcripts in the network.
This means that the large majority of our $(length(vertexlist)) network transcripts are not associated with the significant pathways.
If we can observe how these transcripts interact with the known pathway transcripts, then we can infer functional associations for the novel or less understood transcripts.


	"""
	
end

# ╔═╡ f24286b1-a29f-49dc-8005-058dfcf4440f
begin

	orbit_sigs_file = "$cwd/output/cache/$(params.test_name)/cutoff/$(params.expression_cutoff)/normalisation/$(params.norm_method)/sampling/$(params.variance_percent)/similarity/$(params.coexpression)/threshold/$(params.threshold)/threshold_method/$(params.threshold_method)/analysis/graphlets/validation/coincidents/orbit_sigs.jld2"
	orbit_sigs =  cache_load(orbit_sigs_file,"orbit_sigs")
	                        
	                
end;   

# ╔═╡ c5f26162-6c29-46c2-b08b-8832ea301207
begin
	  ## now compare the significance profile of those nodes that are not attached to a pathway to the average pathway profile of known pathway nodes
		                ##convert to array form for comparisons
		                orbit_sigs_array = map(x->Array(x[!,2:end]),orbit_sigs)

end;

# ╔═╡ 33367e44-3233-4200-a091-a4bd4ea9adeb
begin


	output_dir = "$cwd/output/plots/pluto_test"
	run(`mkdir -p $output_dir`)
	                # Gadfly beeswarm visualisation:
	        		# get data into wide format
	                # size of each df
	                last_col = size(orbit_sigs[1])[2]-1
	                wide_orbit_sigs = vcat(map(x->stack(x,2:last_col+1),orbit_sigs)...)
	                #for each pathway, we map the three orbit categories side by side
	                palette = ["#db63c5","#bababa","#32a852"]
					plot_grid = Array{Plot,1}(undef,length(candidate_pathways)) 
	    
		for (i,p) in enumerate(candidate_pathways)
	                    @info "Drawing beeswarm for $p..."    
	                    p_df = filter(:Pathway=>x->x == p,wide_orbit_sigs)
					
	                    #insertcols!(p_df,:log_value =>log.(p_df.value))
	                    #define colors by whether entrez id of node is in pathway
	                    in_pathway = vcat(collect(eachrow(repeat(in.(1:length(vertexlist),Ref(candidates[p])),1,last_col)))...)
	                    non_entrez = vcat(collect(eachrow(repeat(entrez_id_vector.==0,1,last_col)))...)
	                    coloring = CategoricalArray(in_pathway - non_entrez)
	                  coloring = recode(coloring,-1=>"unidentified",0=>"not in pathway",1=>"in pathway")
	                    insertcols!(p_df,:color =>coloring)
	                    #remove non pathway nodes
	                    filter!(:color=>x->x=="in pathway",p_df)
	                    #remove zero nodes
	                    #filter!(:value=>x->x!=0,p_df)
			
	                    pl = plot(p_df,x = :variable,y = :value, color = :color,Guide.title(p),Geom.beeswarm(padding = 1mm),Theme(bar_spacing=1mm,point_size=0.5mm),Scale.color_discrete_manual(palette...));
	                    plot_grid[i] = pl
						#draw(SVG("$(output_dir)/$(p)_beeswarm.svg",30cm,20cm),pl)
	                end
	
end;

# ╔═╡ 812abebd-a67d-498c-939b-98f9b4d8a300
plot_grid[findall(x-> occursin("Glyoxylate", x),candidate_pathways)...]

# ╔═╡ b0ce6943-a214-4f5c-a75a-c3d59c71aa81
plot_grid[findall(x-> occursin("Alzheimer", x),candidate_pathways)...]

# ╔═╡ 92c7926c-f21f-4665-84df-1ea98229cd4d
begin 	
	                #High zero exploration
                ## over all transcripts 
	                total_zero_proportion = sum(map(x->x.==0,orbit_sigs_array))./length(orbit_sigs_array)
	            
	                
	# #compare for all pathways,speficically for that pathway
	                zero_scores = zeros(Int,length(candidate_pathways),last_col) 
	                for (i,p) in enumerate(candidate_pathways)
	                    subset = candidates[p] 
	                    subset_zero_proportion = sum(map(x->x.==0,orbit_sigs_array[subset]))./length(orbit_sigs_array[subset])
	                #     #compare
	                    zero_scores[i,:] = (subset_zero_proportion.<total_zero_proportion)[i,:]
	                end
	                # #find those pathways with majority of zero proportions below total proportions 
	                zero_passes = vec(sum(zero_scores,dims=2).>(last_col/2))
	                zero_candidate_pathways = candidate_pathways[zero_passes]
	                zero_orbit_sigs = map(x->filter(:Pathway=>y->y in zero_candidate_pathways,x),orbit_sigs)
	                zero_orbit_sigs_array = map(x->Array(x[!,2:end]),zero_orbit_sigs)
	                zero_candidates = Dict(Pair.(zero_candidate_pathways,[candidates[x] for x in zero_candidate_pathways]))
  
end;

# ╔═╡ f7b11d20-9e7e-4197-b9b4-136d3e86644f
DataFrame(pathways = zero_candidate_pathways)

# ╔═╡ df210827-7a7e-4181-bdae-c65974ff81db
zero_candidate_pathways[8]

# ╔═╡ 534b73a8-ced9-4e8c-a9df-ef73af97198d
begin
                #ecdfs
                #store ecdf functions in table
                ecdf_table = Array{ECDF,2}(undef,length(zero_candidate_pathways),last_col)
                for (i,p) in enumerate(zero_candidate_pathways)
                    for j in 1:last_col
                        ecdf_table[i,j] = ecdf(map(x->x[i,j],zero_orbit_sigs_array))
                    end
                end
                #first check candidate probabilities across all categories
                known_pathway_probs = Array{Array{Float64,2},1}(undef,length(zero_candidate_pathways))
                for (i,c) in enumerate(zero_candidate_pathways)
                    known_pathway_probs[i] = hcat([map(ecdf_table[i,j],map(x->x[i,j],orbit_sigs_array[zero_candidates[c]])) for j in 1:last_col]...)
                end

end;

# ╔═╡ 87c1ea4a-26b2-42af-a21d-b6ff80366562
begin
	                sig_pathway_occurences = countmap(vcat([candidates[x] for x in zero_candidate_pathways]...))

	orbit_names = names(orbit_sigs[1])[2:last_col+1]
                #collect known pathway vectors for corresponding known pathway nodes
                known_pathway_dfs = Array{DataFrame,1}(undef,length(zero_candidate_pathways))
                for (i,p) in enumerate(zero_candidate_pathways)
                    subset = zero_candidates[p]
                    df_build = DataFrame(hcat(map(x->x[i,:],zero_orbit_sigs_array[subset])...)',orbit_names)
                    insertcols!(df_build,1,:shared=>[sig_pathway_occurences[x] for x in subset].-1)
                    insertcols!(df_build,1,:transcript_id=>subset)
                    known_pathway_dfs[i] = df_build 
                    #print("Known nodes in pathway $p are in this many other pathways:\n")
                    #print([sig_pathway_occurences[x] for x in subset].-1)
                    #print("\n")
                end
                known_pathway_arrays = map(x->Array(x[!,3:end]),known_pathway_dfs)
                #known ecdfs
                #calculate ecdfs just for known pathway node values (each row corresponds to a pathway, each column to an orbit) 
                known_ecdf_table = Array{ECDF,2}(undef,length(zero_candidate_pathways),last_col)
                for (i,p) in enumerate(zero_candidate_pathways)
                    known_ecdf_table[i,:] = ecdf.(eachcol(known_pathway_arrays[i]))
                end
end


  

# ╔═╡ c2d2fe61-9e00-4db8-aa80-97deb551486a
begin
	                #shape plots into a grid
                ncols = 6
                dims = fldmod(length(zero_candidate_pathways),ncols)
                ecdf_plots = Array{Union{Plot,Context},2}(undef,dims[1]+(dims[2]>0),ncols)
                #plots = Array{Union{Plot,Context},1}(undef,length(keys(candidates)))
                for (i,p) in enumerate(zero_candidate_pathways)
                    #find max over all measures for pathway
                               m = max(known_pathway_arrays[i]...)
                    ecdf_plots[i] = plot([layer(x->known_ecdf_table[i,j](x),0,m,color=[orbit_names[j]]) for j in 1:last_col]...,
                              Scale.color_discrete_manual("orange", "green", "purple"),
                              Guide.title(p)
                              ,Guide.xlabel("count")
                              #,Theme(major_label_font_size=4pt,key_position=:none)
						    ,Guide.colorkey(title="orbit position")
					);
                      
                              #Guide.title(p));
                end
               #Add legend pane
               legend = plot(wide_orbit_sigs,color=:variable,
                             Geom.blank,
                             Scale.color_discrete_manual("orange", "green", "purple"));
                #append blank gridspots if necessary
                for i in 1:length(ecdf_plots)
                    if(!isassigned(ecdf_plots,i))
                        ecdf_plots[i] = context()
                    end
                end
                #TODO this needs to be more generalised for any dimension of gridstack/number of pathways 
                ecdf_plots[2,6] = legend;
                #draw(SVG("$(output_dir)/known_ecdfs.svg",30cm,20cm),gridstack(plots))

						end;

# ╔═╡ 18dbd2ba-ad4c-40fc-92bc-3b8579d8d952
ecdf_plots[1,2]

# ╔═╡ 2c0de2b0-3d95-4abd-a7fe-7bbcf3fcfd52
ecdf_plots[1,4]

# ╔═╡ 51ee4885-582a-4349-8d91-98f583b8a510
ecdf_plots[2,4]

# ╔═╡ 3d0da8ba-1b04-462a-82e3-8b50eb1c29d8
begin
	              ##compare each unknown node to known ecdfs
                #set threshold for which to check on i.e. is node orbit count higher than the top ~thresh~% of known nodes
                thresh = 0.05
                ub = 1 -thresh
                #first check if there are any orbits for any pathways that have such low representation that 0 count exceed threshold. we will filter these out of analysis
                low_filter = .!(map(x->x(0),known_ecdf_table).>ub)
                #for each node, map each known ecdf to the corresponding orbit count and check against threshold, factoring in the low filter
                unknown_ecdf_comparison = map(y->(reshape(map(x->known_ecdf_table[x](y[x]),1:length(known_ecdf_table)),size(known_ecdf_table)[1],size(known_ecdf_table)[2]).>ub).*low_filter,zero_orbit_sigs_array)
                #determine a node as significantly linked to a pathway if at least half its orbit counts are significant
                sig_check = map(x->(sum(x,dims=2).>(last_col/2)),unknown_ecdf_comparison)
                sig_nodes= findall(x->sum(x)>0,sig_check)
                sig_nodes_dict = Dict(Pair.(sig_nodes,map(x->zero_candidate_pathways[vec(x)],sig_check[sig_nodes])))


end;

# ╔═╡ 7a0255c6-82c6-48c6-83e1-cd83606b05dc
begin
	lf = DataFrame(Array(low_filter),orbit_names)
	
	insertcols!(lf,1,:pathway => zero_candidate_pathways)
end

# ╔═╡ 84aa8941-8766-4482-ba97-cfb5c477a072
DataFrame(node = collect(keys(sig_nodes_dict)), pathways = collect(values(sig_nodes_dict)),type=vertexlist[sig_nodes])

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CategoricalArrays = "324d7699-5711-5eae-9e2f-1d82baa6b597"
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
Compose = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DataPreprocessing = "0c67aaa8-d5ff-4929-99a0-75b09377fbc9"
Gadfly = "c91e804a-d5a3-530f-b6f0-dfbca275c004"
GraphletAnalysis = "32f39a16-8143-4a50-a7e7-080c0e917f42"
GraphletCounting = "7ac45bc0-02f1-46da-ad35-65e91b15b4e1"
JLD = "4138dd39-2aa7-5051-a626-17a0bb65d9c8"
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
NetworkConstruction = "6c2e41d2-72ae-425a-84e9-b8f08a301efb"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
ProjectFunctions = "a8586eae-54f0-4952-9436-ba92c8ab3181"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
CSV = "~0.9.11"
CategoricalArrays = "~0.10.2"
Colors = "~0.12.8"
Compose = "~0.9.2"
DataFrames = "~1.3.0"
DataPreprocessing = "~0.1.0"
Gadfly = "~1.3.4"
GraphletAnalysis = "~0.1.0"
GraphletCounting = "~0.1.0"
JLD = "~0.12.3"
JLD2 = "~0.4.15"
NetworkConstruction = "~0.1.0"
ProjectFunctions = "~0.1.0"
Revise = "~3.1.20"
StatsBase = "~0.33.13"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "485ee0867925449198280d4af84bdb46a2a404d0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.0.1"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
git-tree-sha1 = "bdf73eec6a88885256f282d48eafcad25d7de494"
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "f87e559f87a45bece9c9ed97458d3afe98b1ebb9"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.1.0"

[[Artifacts]]
deps = ["Pkg"]
git-tree-sha1 = "c30985d8821e0cd73870b17b0ed0ce6dc44cb744"
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.3.0"

[[AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Blosc]]
deps = ["Blosc_jll"]
git-tree-sha1 = "217da19d6f3a94753e580a8bc241c7cbefd9281f"
uuid = "a74b3585-a348-5f62-a45c-50e91977d574"
version = "0.7.1"

[[Blosc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Lz4_jll", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "e747dac84f39c62aff6956651ec359686490134e"
uuid = "0b7ba130-8d10-5ba8-a3d6-c5182647fed9"
version = "1.21.0+0"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c3598e525718abcc440f69cc6d5f60dda0a1b61e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.6+5"

[[CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "49f14b6c56a2da47608fe30aed711b5882264d7a"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.9.11"

[[Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "e2f47f6d8337369411569fd45ae5753ca10394c6"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.0+6"

[[CategoricalArrays]]
deps = ["DataAPI", "Future", "Missings", "Printf", "Requires", "Statistics", "Unicode"]
git-tree-sha1 = "c308f209870fdbd84cb20332b6dfaf14bf3387f8"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.10.2"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4c26b4e9e91ca528ea212927326ece5918a04b47"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.2"

[[ChangesOfVariables]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "9a1d594397670492219635b35a3d830b04730d62"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.1"

[[CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "9aa8a5ebb6b5bf469a7e0e2b5202cf6f8c291104"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.0.6"

[[CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dce3e3fea680869eaa0b774b2e8343e9ff442313"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.40.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "8e695f735fca77e9708e795eda62afdb869cbb70"
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.3.4+0"

[[Compose]]
deps = ["Base64", "Colors", "DataStructures", "Dates", "IterTools", "JSON", "LinearAlgebra", "Measures", "Printf", "Random", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "c6461fc7c35a4bb8d00905df7adafcff1fe3a6bc"
uuid = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
version = "0.9.2"

[[Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "6cdc8832ba11c7695f494c9d9a1c31e90959ce0f"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.6.0"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[CoupledFields]]
deps = ["LinearAlgebra", "Statistics", "StatsBase"]
git-tree-sha1 = "6c9671364c68c1158ac2524ac881536195b7e7bc"
uuid = "7ad07ef1-bdf2-5661-9d2b-286fd4296dac"
version = "0.2.0"

[[Crayons]]
git-tree-sha1 = "3f71217b538d7aaee0b69ab47d9b7724ca8afa0d"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.0.4"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "2e993336a3f68216be91eb8ee4625ebbaba19147"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.0"

[[DataPreprocessing]]
deps = ["CSV", "Cairo", "Compose", "DataFrames", "Gadfly", "LinearAlgebra", "RCall", "Statistics"]
path = "../../home/osbornejr/app/packages/DataPreprocessing"
uuid = "0c67aaa8-d5ff-4929-99a0-75b09377fbc9"
version = "0.1.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[Discretizers]]
deps = ["DataStructures", "SpecialFunctions", "Statistics", "StatsBase"]
git-tree-sha1 = "88a41c96d120c5b76a549601e5800c5743715184"
uuid = "6e83dbb3-75ca-525b-8ae2-3751f0dd50b4"
version = "3.2.2"

[[Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "7f3bec11f4bcd01bc1f507ebce5eadf1b0a78f47"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.34"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
git-tree-sha1 = "cd7f38a6540ce18f53dba40f8e251eeab49871c6"
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.5.2"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "1402e52fcda25064f51c77a9655ce8680b76acf0"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.7+6"

[[FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "IntelOpenMP_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Reexport"]
git-tree-sha1 = "1b48dbde42f307e48685fa9213d8b9f8c0d87594"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.3.2"

[[FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3676abafff7e4ff07bbd2c42b3d8201f31653dcc"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.9+8"

[[FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "2db648b6712831ecb333eae76dbfd1c156ca13bb"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.11.2"

[[FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "04d13bfa8ef11720c24e4d840c0033d145537df7"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.17"

[[FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "8756f9935b7ccc9064c6eef0bff0ad643df733a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.7"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "35895cf184ceaab11fd778b4590144034a167a2f"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.1+14"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "cbd58c9deb1d304f5a245a0b7eb841a2560cfec6"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.1+5"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0d20aed5b14dd4c9a2453c1b601d08e1149679cc"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.5+6"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[Gadfly]]
deps = ["Base64", "CategoricalArrays", "Colors", "Compose", "Contour", "CoupledFields", "DataAPI", "DataStructures", "Dates", "Distributions", "DocStringExtensions", "Hexagons", "IndirectArrays", "IterTools", "JSON", "Juno", "KernelDensity", "LinearAlgebra", "Loess", "Measures", "Printf", "REPL", "Random", "Requires", "Showoff", "Statistics"]
git-tree-sha1 = "13b402ae74c0558a83c02daa2f3314ddb2d515d3"
uuid = "c91e804a-d5a3-530f-b6f0-dfbca275c004"
version = "1.3.4"

[[Gettext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "8c14294a079216000a0bdca5ec5a447f073ddc9d"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.20.1+7"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "04690cc5008b38ecbdfede949220bc7d9ba26397"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.59.0+4"

[[Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "1c5a84319923bea76fa145d49e93aa4394c73fc2"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.1"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42adbc6fd39ba41138f894b8ac711146a2b0d986"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.13+4"

[[GraphletAnalysis]]
deps = ["DataFrames", "DataPreprocessing", "GraphletCounting", "LightGraphs", "NetworkConstruction", "RCall"]
path = "../../home/osbornejr/app/packages/GraphletAnalysis"
uuid = "32f39a16-8143-4a50-a7e7-080c0e917f42"
version = "0.1.0"

[[GraphletCounting]]
deps = ["DataFrames", "DataStructures", "Distributed", "LightGraphs", "ProgressMeter", "RCall"]
path = "../../home/osbornejr/app/packages/GraphletCounting"
uuid = "7ac45bc0-02f1-46da-ad35-65e91b15b4e1"
version = "0.1.0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HDF5]]
deps = ["Blosc", "Compat", "HDF5_jll", "Libdl", "Mmap", "Random", "Requires"]
git-tree-sha1 = "698c099c6613d7b7f151832868728f426abe698b"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.15.7"

[[HDF5_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "OpenSSL_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "fd83fa0bde42e01952757f01149dd968c06c4dba"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.12.0+1"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Gettext_jll", "Glib_jll", "Graphite2_jll", "ICU_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "90bed5fc61d12d10832ebf988988104888eebaca"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.6.1+10"

[[Hexagons]]
deps = ["Test"]
git-tree-sha1 = "de4a6f9e7c4710ced6838ca906f81905f7385fd6"
uuid = "a1b4810d-1bce-5fbd-ac56-80944d57a21f"
version = "0.2.0"

[[ICU_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ce9cf3c5490b045710d60ac3fd2fe48188846b3"
uuid = "a51ab1cf-af8e-5615-a023-bc2c838bba6b"
version = "67.1.0+3"

[[IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[InformationMeasures]]
deps = ["Discretizers"]
git-tree-sha1 = "874d48f2026e8faf3fd55c86973fd028b02cd1a0"
uuid = "96684042-fbdc-5399-9b8e-d34e539a126c"
version = "0.3.1"

[[InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "ca99cac337f8e0561c6a6edeeae5bf6966a78d21"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.0"

[[IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "61aa005707ea2cebf47c8d780da8dc9bc4e0c512"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.4"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLD]]
deps = ["FileIO", "HDF5", "Printf"]
git-tree-sha1 = "1d291ba1730de859903b480e6f85a0dc40c19dcb"
uuid = "4138dd39-2aa7-5051-a626-17a0bb65d9c8"
version = "0.12.3"

[[JLD2]]
deps = ["DataStructures", "FileIO", "MacroTools", "Mmap", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "46b7834ec8165c541b0b5d1c8ba63ec940723ffb"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.15"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "e273807f38074f033d94207a201e6e827d8417db"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.8.21"

[[Juno]]
deps = ["Base64", "Logging", "Media", "Profile"]
git-tree-sha1 = "07cb43290a840908a771552911a6274bc6c072c7"
uuid = "e5e0dc1b-0480-54bc-9374-aad01c23163d"
version = "0.8.4"

[[KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "591e8dc09ad18386189610acafb970032c519707"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.3"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f128cd6cd05ffd6d3df0523ed99b90ff6f9b349a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.0+3"

[[LazyArtifacts]]
deps = ["Pkg"]
git-tree-sha1 = "4bb5499a1fc437342ea9ab7e319ede5a457c0968"
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.3.0"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
git-tree-sha1 = "cdbe7465ab7b52358804713a53c7fe1dac3f8a3f"
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[LibCURL_jll]]
deps = ["LibSSH2_jll", "Libdl", "MbedTLS_jll", "Pkg", "Zlib_jll", "nghttp2_jll"]
git-tree-sha1 = "897d962c20031e6012bba7b3dcb7a667170dad17"
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.70.0+2"

[[LibGit2]]
deps = ["Printf"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Libdl", "MbedTLS_jll", "Pkg"]
git-tree-sha1 = "717705533148132e5466f2924b9a3657b16158e8"
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.9.0+3"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "a2cd088a88c0d37eef7d209fd3d8712febce0d90"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.1+4"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "b391a18ab1170a2e568f9fb8d83bc7c780cb9999"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.5+4"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ec7f2e8ad5c9fa99fc773376cdbc86d9a5a23cb7"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.36.0+3"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cba7b560fcc00f8cd770fa85a498cbc1d63ff618"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.0+8"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51ad0c01c94c1ce48d5cad629425035ad030bfd5"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.34.0+3"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f879ae9edbaa2c74c922e8b85bb83cc84ea1450b"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.34.0+7"

[[LightGraphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "432428df5f360964040ed60418dd5601ecd240b6"
uuid = "093fc24a-ae57-5d10-9952-331d41423f4d"
version = "1.3.5"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Loess]]
deps = ["Distances", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "46efcea75c890e5d820e670516dc156689851722"
uuid = "4345ca2d-374a-55d4-8d30-97f9976e7612"
version = "0.5.4"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "be9eef9f9d78cecb6f262f3c10da151a6c5ab827"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.5"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "491a883c4fef1103077a7f648961adbf9c8dd933"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.1.2"

[[Lz4_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "51b1db0732bbdcfabb60e36095cc3ed9c0016932"
uuid = "5ced341a-0733-55b8-9ab6-a4889d929147"
version = "1.9.2+2"

[[MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "5455aef09b40e5020e1520f551fa3135040d4ed0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2021.1.1+2"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0eef589dd1c26a3ac9d753fe1a8bcad63f956fa6"
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.16.8+1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Media]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "75a54abd10709c01f1b86b84ec225d26e840ed58"
uuid = "e89f7d12-3494-54d1-8411-f7d8b9ae1f27"
version = "0.5.0"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f1662575f7bf53c73c2bbc763bace4b024de822c"
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2021.1.19+0"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkConstruction]]
deps = ["DataFrames", "DataPreprocessing", "Distributed", "InformationMeasures", "LightGraphs", "LinearAlgebra", "Printf", "ProgressMeter", "RCall", "SharedArrays", "Statistics", "StatsBase"]
path = "../../home/osbornejr/app/packages/NetworkConstruction"
uuid = "6c2e41d2-72ae-425a-84e9-b8f08a301efb"
version = "0.1.0"

[[NetworkOptions]]
git-tree-sha1 = "ed3157f48a05543cce9b241e1f2815f7e843d96e"
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[OpenLibm_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "d22054f66695fe580009c09e765175cbf7f13031"
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.7.1+0"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "71bbbc616a1d710879f5a1021bcba65ffba6ce58"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.1+6"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9db77584158d0ab52307f8c04f8e7c08ca76b5b3"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.3+4"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "1b556ad51dceefdbf30e86ffa8f528b73c7df2bb"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.42.0+4"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "ee26b350276c51697c9c2d88a072b339f9f03d73"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.5"

[[Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9a336dee51d20d1ed890c4a8dca636e86e2b76ca"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.42.4+10"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "ae4bbcadb2906ccc085cf52ac286dc1377dceccc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.2"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6a20a83c1ae86416f0a5de605eaea08a552844a3"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.0+0"

[[Pkg]]
deps = ["Dates", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "UUIDs"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "db3a23166af8aebf4db5ef87ac5b00d36eb771e2"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.0"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "d940010be611ee9d67064fe559edbb305f8cc0eb"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.2.3"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "afadeba63d90ff223a6a48d2009434ecee2ec9e8"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.1"

[[ProjectFunctions]]
deps = ["JLD2"]
path = "../../home/osbornejr/app/packages/ProjectFunctions"
uuid = "a8586eae-54f0-4952-9436-ba92c8ab3181"
version = "0.1.0"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[RCall]]
deps = ["CategoricalArrays", "Conda", "DataFrames", "DataStructures", "Dates", "Libdl", "Missings", "REPL", "Random", "Requires", "StatsModels", "WinReg"]
git-tree-sha1 = "80a056277142a340e646beea0e213f9aecb99caa"
uuid = "6f49c342-dc21-5d91-9882-a32aef131414"
version = "0.13.12"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Ratios]]
deps = ["Requires"]
git-tree-sha1 = "01d341f502250e81f6fec0afe662aa861392a3aa"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.2"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "41deb3df28ecf75307b6e492a738821b031f8425"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.1.20"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "86c5647b565873641538d8f812c04e4c9dbeb370"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.6.1"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "1b7bf41258f6c5c9c31df8c1ba34c1fc88674957"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.2.2+2"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "f45b34656397a1f6e729901dc9ef679610bd12b5"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.8"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[ShiftedArrays]]
git-tree-sha1 = "22395afdcf37d6709a5a0766cc4a5ca52cb85ea0"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "1.0.0"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f0bccf98e16759818ffc5d97ac3ebf87eb950150"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.8.1"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "0f2aa8e32d511f758a2ce49208181f7733a0936a"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.1.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "2bb0cb32026a66037360606510fca5984ccc6b75"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.13"

[[StatsFuns]]
deps = ["ChainRulesCore", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "bedb3e17cc1d94ce0e6e66d3afa47157978ba404"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.14"

[[StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "REPL", "ShiftedArrays", "SparseArrays", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "677488c295051568b0b79a77a8c44aa86e78b359"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.6.28"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[TOML]]
deps = ["Dates"]
git-tree-sha1 = "44aaac2d2aec4a850302f9aa69127c74f0c3787e"
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "fed34d0e71b91734bf0a7e10eb1bb05296ddbcd0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.0"

[[Test]]
deps = ["Distributed", "InteractiveUtils", "Logging", "Random"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[VersionParsing]]
git-tree-sha1 = "e575cf85535c7c3292b4d89d89cc29e8c3098e47"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.2.1"

[[WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "c69f9da3ff2f4f02e811c3323c22e5dfcb584cfa"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.1"

[[WinReg]]
deps = ["Test"]
git-tree-sha1 = "808380e0a0483e134081cc54150be4177959b5f4"
uuid = "1b915085-20d7-51cf-bf83-8f477d6f5128"
version = "0.3.1"

[[WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "be0db24f70aae7e2b89f2f3092e93b8606d659a6"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.10+3"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "2b3eac39df218762d2d005702d601cd44c997497"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.33+4"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "320228915c8debb12cb434c59057290f0834dbf6"
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.11+18"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "2c1332c54931e83f8f94d310fa447fd743e8d600"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.4.8+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "6abbc424248097d69c0c87ba50fcb0753f93e0ee"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.37+6"

[[nghttp2_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "8e2c44ab4d49ad9518f359ed8b62f83ba8beede4"
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.40.0+2"
"""

# ╔═╡ Cell order:
# ╟─74c92bd0-ef4f-41a0-bc4b-5063ffcd12df
# ╟─940181d4-a9b0-47e4-a13d-db2eb175e22e
# ╟─72b2183e-99a9-4c9b-b74b-7e7966eb4bb8
# ╟─28836056-6604-4918-9f74-39bf81ad0559
# ╟─8267a9ad-9551-4842-9194-5e250e79305e
# ╟─2733aa29-9e63-4b60-a7ee-0e9abcd62974
# ╟─e01833de-63e0-493a-8d6c-0240c1fe1633
# ╟─3c79ec7e-7e64-4f24-aa65-35c79167301b
# ╠═08edac1b-f869-4464-828d-bcb60c64a98e
# ╟─752f5e2a-6a63-442f-9e28-90db65b6eeb1
# ╟─f759d399-f76a-4cb5-b0d6-bf64aab1380a
# ╟─a380f1ed-32bb-46f1-bd2c-8c633afa70d0
# ╟─6a2d3d61-cc8e-48ba-8839-705979033484
# ╟─812abebd-a67d-498c-939b-98f9b4d8a300
# ╟─b0ce6943-a214-4f5c-a75a-c3d59c71aa81
# ╟─2abc718c-e711-4b7f-880c-38922e225dd2
# ╟─4ee8ef83-ee2a-458c-b7ce-bf16fd4c5baa
# ╟─f7b11d20-9e7e-4197-b9b4-136d3e86644f
# ╟─822164ee-0193-4c16-a245-b55ea8082530
# ╟─18dbd2ba-ad4c-40fc-92bc-3b8579d8d952
# ╟─2c0de2b0-3d95-4abd-a7fe-7bbcf3fcfd52
# ╟─51ee4885-582a-4349-8d91-98f583b8a510
# ╟─df210827-7a7e-4181-bdae-c65974ff81db
# ╟─07904aca-d5da-444a-8fa5-55f6af49eb23
# ╟─7a0255c6-82c6-48c6-83e1-cd83606b05dc
# ╟─7f35a6f0-71e3-4448-8dd1-fe359ac73f94
# ╠═84aa8941-8766-4482-ba97-cfb5c477a072
# ╟─877da2cf-fe40-4bfd-b44c-478427f3c81f
# ╟─99f1209a-3dea-4c60-8b5b-fe7da9219a1d
# ╟─d2f8505a-6d00-4a75-85d9-aba9bb251302
# ╟─3f3e6d45-d576-4385-bea8-e55a37d34512
# ╟─f34235b1-0e47-4f28-94db-ce3cfa598a91
# ╟─9bb259ec-8528-4049-80bf-5fa0d543e47c
# ╟─f24286b1-a29f-49dc-8005-058dfcf4440f
# ╟─c5f26162-6c29-46c2-b08b-8832ea301207
# ╟─33367e44-3233-4200-a091-a4bd4ea9adeb
# ╟─92c7926c-f21f-4665-84df-1ea98229cd4d
# ╟─534b73a8-ced9-4e8c-a9df-ef73af97198d
# ╟─c2d2fe61-9e00-4db8-aa80-97deb551486a
# ╟─87c1ea4a-26b2-42af-a21d-b6ff80366562
# ╟─3d0da8ba-1b04-462a-82e3-8b50eb1c29d8
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
