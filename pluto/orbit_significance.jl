### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 74c92bd0-ef4f-41a0-bc4b-5063ffcd12df
# ╠═╡ show_logs = false
### set up notebook environment
begin
cwd = ENV["PWD"];	
	import Pkg
	dev = Pkg.develop
	dev(path=cwd*"/packages/DataPreprocessing")
	dev(path=cwd*"/packages/NetworkConstruction")
	dev(path=cwd*"/packages/GraphletCounting")
	dev(path=cwd*"/packages/GraphletAnalysis")
	dev(path=cwd*"/packages/ProjectFunctions")
	using Revise	
	using CommonMark
	using PlutoUI
	using JLD2
	using StatsBase
	using Gadfly,CategoricalArrays,Colors, Compose
	using DataFrames
	using YAML
	using Infiltrator
	## load dev packages last
	using DataPreprocessing
	using NetworkConstruction
	using GraphletCounting
	using GraphletAnalysis
	using ProjectFunctions
	
	TableOfContents()
end;

# ╔═╡ 6ea01822-2aa0-4073-84b0-304e5a86f9ea
TableOfContents()

# ╔═╡ 940181d4-a9b0-47e4-a13d-db2eb175e22e
# ╠═╡ show_logs = false
### Setup input data (from human smoker GSE68559 dataset) loaded from existing cache
begin
	#using Logging
	#redirect_stdout(devnull) do
		#logger = ConsoleLogger(stdout)
		#with_logger(logger) do
			config_file = "$cwd/config/run-files/pluto.yaml"
			ProjectFunctions.load_config(config_file)
			##move params values to specific variables for now (easier)
			params_test_name = ProjectFunctions.params["test_name"]
			params_null_model_size = ProjectFunctions.params["analysis"]["null_model_size"] 
			params_expression_cutoff = ProjectFunctions.params["data_preprocessing"]["expression_cutoff"] 
			params_variance_percent = ProjectFunctions.params["data_preprocessing"]["variance_percent"] 
			params_norm_method = ProjectFunctions.params["data_preprocessing"]["norm_method"] 
			params_threshold = ProjectFunctions.params["network_construction"]["threshold"] 
			params_threshold_method = ProjectFunctions.params["network_construction"]["threshold_method"] 
			params_coexpression = ProjectFunctions.params["network_construction"]["coexpression"]
			params_graphlet_size = ProjectFunctions.params["analysis"]["graphlet_size"]
			raw_counts = ProjectFunctions.get_input_data();
			processed_counts = ProjectFunctions.data_preprocessing(raw_counts);
			    adj_matrix,network_counts,vertexlist,edgelist = ProjectFunctions.network_construction(processed_counts)
		
			#biomart_raw_counts_file = "$cwd/output/cache/$(params_test_name)_raw_counts_biomart.jld2";
			#raw_counts = ProjectFunctions.cache_load(biomart_raw_counts_file,"raw counts");
		#end
	#end
	
end;

# ╔═╡ 3f3e6d45-d576-4385-bea8-e55a37d34512
begin
	using CSV
	coincidents_file = "$cwd/output/cache/$(params_test_name)/expression_cutoff/$(params_expression_cutoff)/normalisation/$(params_norm_method)/sampling/$(params_variance_percent)/similarity/$(params_coexpression)/threshold/$(params_threshold)/threshold_method/$(params_threshold_method)/analysis/graphlets/$(params_graphlet_size)/coincidents/coincidents.csv";
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

# ╔═╡ 312f20e3-3c97-4108-847d-f2276abfab9a
begin
	degrees = vec(sum(adj_matrix,dims=2))
	#degrees = sort(vec(degrees))
	p = plot(x = 1:length(degrees),y=degrees, Geom.bar,Guide.xlabel("id"),Guide.ylabel("count"))

end;

# ╔═╡ 7211b91e-3ca1-4688-b88f-e9e5e92ca918
Coincidents

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
"""

# ╔═╡ 67f24598-4b27-4dd8-b533-45c1cc4e44a6
md"""
## Network details

The network we are working with is derived from a *homo sapiens* dataset[^1] that compared 10 different regions of brain tissue across smokers and nonsmokers.
Transcripts are classified as either `coding` or `noncoding` based on a multistep process that factors in coding potential and checks for preexisting records on the Ensembl database.
After normalisation, transcripts are selected as potential nodes if they show enough variable expression across samples.
Edges between nodes are then defined using a partial information measure method that focuses on identify only the **direct** interactions between transcripts.[^2]
The resulting network has the following structure:
"""



# ╔═╡ 4861b8cf-cf59-4cfd-b412-089ccfc00e90
md"""
|nodes | noncoding nodes | coding nodes | number of edges |
|:----:|:---------------:|:------------:|:---------------:|
|$(length(vertexlist))   | $(length(findall(x->x=="noncoding",vertexlist)))  |$(length(findall(x->x=="coding",vertexlist)))   |$(length(edgelist))   |      

"""

# ╔═╡ df21e8f4-7af6-4401-a8d7-e4f21e9c7139
begin
	parser = enable!(Parser(),[FootnoteRule()])
	
end;

# ╔═╡ 4ef32b7a-d7f3-4642-a75c-2870a5130ade
md"""
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

# ╔═╡ 0c721876-736e-495e-85b6-490ea7ded8f2


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

# ╔═╡ f34235b1-0e47-4f28-94db-ce3cfa598a91
begin
	
		
	    sub_Coincidents = filter(:Hom_graphlet=>x->occursin("4-",x),Coincidents);
	
end;

# ╔═╡ 08edac1b-f869-4464-828d-bcb60c64a98e
sub_Coincidents

# ╔═╡ 9bb259ec-8528-4049-80bf-5fa0d543e47c
begin
	kegg_file = "/home/osbornejr/app/output/cache/GSE68559/expression_cutoff/$(params_expression_cutoff)/normalisation/$(params_norm_method)/sampling/$(params_variance_percent)/similarity/$(params_coexpression)/threshold/$(params_threshold)/threshold_method/$(params_threshold_method)/analysis/graphlets/$(params_graphlet_size)/coincidents/kegg_info.jld2"

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

# ╔═╡ c2ddfc1b-5ef7-4c0f-82fb-d86c015aba74
cm"""
## References
[^1]: @WebbRNAsequencingtranscriptomes2015
[^2]: @ChanGeneRegulatoryNetwork2017
"""

# ╔═╡ f24286b1-a29f-49dc-8005-058dfcf4440f
begin

	orbit_sigs_file = "$cwd/output/cache/$(params_test_name)/expression_cutoff/$(params_expression_cutoff)/normalisation/$(params_norm_method)/sampling/$(params_variance_percent)/similarity/$(params_coexpression)/threshold/$(params_threshold)/threshold_method/$(params_threshold_method)/analysis/graphlets/$(params_graphlet_size)/coincidents/orbit_sigs.jld2"
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
			
	                    pl = plot(p_df,x = :variable,y = :value, color = :color,Guide.title(p),Guide.xlabel("Orbit"),Guide.ylabel("Orbit frequency (per node)",orientation=:vertical),Geom.beeswarm(padding = 1mm),Theme(bar_spacing=1mm,point_size=0.5mm),Scale.color_discrete_manual(palette...));
	                    plot_grid[i] = pl
						#draw(SVG("$(output_dir)/$(p)_beeswarm.svg",30cm,20cm),pl)
	                end
	
end;

# ╔═╡ 812abebd-a67d-498c-939b-98f9b4d8a300
begin
	plot_1 = "Thermogenesis" 
	plot_grid[findall(x-> occursin(plot_1, x),candidate_pathways)...]
end

# ╔═╡ b0ce6943-a214-4f5c-a75a-c3d59c71aa81
begin
	plot_2 = "Alzheimer disease" 
	plot_grid[findall(x-> occursin(plot_2, x),candidate_pathways)...]
end

# ╔═╡ 2abc718c-e711-4b7f-880c-38922e225dd2
md"""
Each node  that is known to be in the pathway is represented as a point for each orbit.
The height of the point indicates the number of times that node occurred in that orbit.
In the first pathway, *$(plot_1)*, there is little  signal across any of the orbits. 
In contrast, the *$(plot_2)* pathway clearly has a lot of coincidence at the graphlet level.
We can deduce that pathways with a high graphlet connectivity in our network are those that are most active as the known pathway nodes are interacting with each other.
This means we can then identify other nodes that are also have a interaction with these known pathway nodes and infer that these nodes also have some importance in the KEGG pathway. 
Thus we narrow our scope to only those pathways that have a high connectivity across orbit categories.
The pathways we exclude may have a significant number of transcripts present in the network, but (at least at the four node level) they do not appear to be active across the sample space. 
"""


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
	                ##quick fix to remove above zero score process
					zero_candidate_pathways = candidate_pathways
	                zero_orbit_sigs = map(x->filter(:Pathway=>y->y in zero_candidate_pathways,x),orbit_sigs)
	                zero_orbit_sigs_array = map(x->Array(x[!,2:end]),zero_orbit_sigs)
	                zero_candidates = Dict(Pair.(zero_candidate_pathways,[candidates[x] for x in zero_candidate_pathways]))
  
end;

# ╔═╡ f7b11d20-9e7e-4197-b9b4-136d3e86644f
DataFrame(pathways = zero_candidate_pathways)

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
ecdf_plots[3,3]

# ╔═╡ 8cdaf845-43c6-4455-822f-d30fbd2ff19b
ecdf_plots[1,3]

# ╔═╡ 3d0da8ba-1b04-462a-82e3-8b50eb1c29d8
begin
	              ##compare each unknown node to known ecdfs
                #set threshold for which to check on i.e. is node orbit count higher than the top ~thresh~% of known nodes
                thresh = 0.05
                ub = 1 -thresh
                #first check if there are any orbits for any pathways that have such low representation that 0 count exceed threshold. we will filter these out of analysis
                low_filter = .!(map(x->x(0),known_ecdf_table).>ub)
				#low_filter = BitMatrix(ones(Bool,20,11))
                #for each node, map each known ecdf to the corresponding orbit count and check against threshold, factoring in the low filter
                unknown_ecdf_comparison = map(y->(reshape(map(x->known_ecdf_table[x](y[x]),1:length(known_ecdf_table)),size(known_ecdf_table)[1],size(known_ecdf_table)[2]).>ub).*low_filter,zero_orbit_sigs_array)
                #determine a node as significantly linked to a pathway if at least half its orbit counts are significant
				sig_cut = 1
                sig_check = map(x->(sum(x,dims=2).>(sig_cut)),unknown_ecdf_comparison)
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

# ╔═╡ ebaa8810-9096-4cc2-9c28-31143143bf14
unique(vcat(values(sig_nodes_dict)...))

# ╔═╡ 61a6f729-141c-4892-b983-dd48d78b6b7d
unique(vcat(values(sig_nodes_dict)...))

# ╔═╡ 726ac732-b9c4-439a-867c-61fcc34abedb


# ╔═╡ f520c106-6e0d-49ec-8b9d-d656daa5b8dd
orbit_sigs_array

# ╔═╡ ea2da746-ac59-4a54-99c1-5ab39a2f3221
length(sig_nodes_dict)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CategoricalArrays = "324d7699-5711-5eae-9e2f-1d82baa6b597"
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
CommonMark = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
Compose = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DataPreprocessing = "0c67aaa8-d5ff-4929-99a0-75b09377fbc9"
Gadfly = "c91e804a-d5a3-530f-b6f0-dfbca275c004"
GraphletAnalysis = "32f39a16-8143-4a50-a7e7-080c0e917f42"
GraphletCounting = "7ac45bc0-02f1-46da-ad35-65e91b15b4e1"
Infiltrator = "5903a43b-9cc3-4c30-8d17-598619ec4e9b"
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
NetworkConstruction = "6c2e41d2-72ae-425a-84e9-b8f08a301efb"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
ProjectFunctions = "a8586eae-54f0-4952-9436-ba92c8ab3181"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
YAML = "ddb6d928-2868-570f-bddf-ab3f9cf99eb6"

[compat]
CSV = "~0.10.9"
CategoricalArrays = "~0.10.7"
Colors = "~0.12.10"
CommonMark = "~0.8.11"
Compose = "~0.9.5"
DataFrames = "~1.5.0"
DataPreprocessing = "~0.1.0"
Gadfly = "~1.3.4"
GraphletAnalysis = "~0.1.0"
GraphletCounting = "~0.1.0"
Infiltrator = "~1.6.3"
JLD2 = "~0.4.31"
NetworkConstruction = "~0.1.0"
PlutoUI = "~0.7.50"
ProjectFunctions = "~0.1.0"
Revise = "~3.5.2"
StatsBase = "~0.33.21"
YAML = "~0.4.8"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "c587dc69b416d33d3eb5703dd395c45afc38a210"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "16b6dbc4cf7caee4e1e75c49485ec67b667098a0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.3.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cc37d689f599e8df4f464b2fa3870ff7db7492ef"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "SnoopPrecompile", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "c700cce799b51c9045473de751e9319bdd1c6e94"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.9"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CategoricalArrays]]
deps = ["DataAPI", "Future", "Missings", "Printf", "Requires", "Statistics", "Unicode"]
git-tree-sha1 = "5084cc1a28976dd1642c9f337b28a3cb03e0f7d2"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.10.7"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "485193efd2176b88e6622a39a246f8c5b600e74e"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.6"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "d730914ef30a06732bdd9f763f6cc32e92ffbff1"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.CommonMark]]
deps = ["Crayons", "JSON", "SnoopPrecompile", "URIs"]
git-tree-sha1 = "4a52799aee66e9528bd59e9c0bdd9322c0bf6998"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.11"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.Compose]]
deps = ["Base64", "Colors", "DataStructures", "Dates", "IterTools", "JSON", "LinearAlgebra", "Measures", "Printf", "Random", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "bf6570a34c850f99407b494757f5d7ad233a7257"
uuid = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
version = "0.9.5"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "e32a90da027ca45d84678b826fffd3110bb3fc90"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.8.0"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.CoupledFields]]
deps = ["LinearAlgebra", "Statistics", "StatsBase"]
git-tree-sha1 = "6c9671364c68c1158ac2524ac881536195b7e7bc"
uuid = "7ad07ef1-bdf2-5661-9d2b-286fd4296dac"
version = "0.2.0"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "aa51303df86f8626a962fccb878430cdb0a97eee"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.5.0"

[[deps.DataPreprocessing]]
deps = ["CSV", "Cairo", "Compose", "DataFrames", "Gadfly", "LinearAlgebra", "RCall", "Statistics"]
path = "/Users/osbornejr/git/heterogeneous-graphlet-counting/packages/DataPreprocessing"
uuid = "0c67aaa8-d5ff-4929-99a0-75b09377fbc9"
version = "0.1.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.Discretizers]]
deps = ["DataStructures", "SpecialFunctions", "Statistics", "StatsBase"]
git-tree-sha1 = "88a41c96d120c5b76a549601e5800c5743715184"
uuid = "6e83dbb3-75ca-525b-8ae2-3751f0dd50b4"
version = "3.2.2"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "49eba9ad9f7ead780bfb7ee319f962c811c6d3b2"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.8"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "13027f188d26206b9e7b863036f87d2f2e7d013a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.87"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "f9818144ce7c8c41edf5c4c179c684d92aa4d9fe"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.6.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "7be5f99f7d15578798f338f5433b6c432ea8037b"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "e27c4ebe80e8699540f2d6c805cc12203b614f12"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.20"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "fc86b4fd3eff76c3ce4f5e96e2fdfa6282722885"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.0.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.Gadfly]]
deps = ["Base64", "CategoricalArrays", "Colors", "Compose", "Contour", "CoupledFields", "DataAPI", "DataStructures", "Dates", "Distributions", "DocStringExtensions", "Hexagons", "IndirectArrays", "IterTools", "JSON", "Juno", "KernelDensity", "LinearAlgebra", "Loess", "Measures", "Printf", "REPL", "Random", "Requires", "Showoff", "Statistics"]
git-tree-sha1 = "13b402ae74c0558a83c02daa2f3314ddb2d515d3"
uuid = "c91e804a-d5a3-530f-b6f0-dfbca275c004"
version = "1.3.4"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.GraphPlot]]
deps = ["ArnoldiMethod", "ColorTypes", "Colors", "Compose", "DelimitedFiles", "Graphs", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "5cd479730a0cb01f880eff119e9803c13f214cab"
uuid = "a2cc645c-3eea-5389-862e-a155d0052231"
version = "0.5.2"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.GraphletAnalysis]]
deps = ["DataFrames", "DataPreprocessing", "GraphletCounting", "Graphs", "NetworkConstruction", "RCall"]
path = "/Users/osbornejr/git/heterogeneous-graphlet-counting/packages/GraphletAnalysis"
uuid = "32f39a16-8143-4a50-a7e7-080c0e917f42"
version = "0.1.0"

[[deps.GraphletCounting]]
deps = ["CSV", "DataFrames", "DataStructures", "Distributed", "Graphs", "LinearAlgebra", "ProgressMeter", "StatsBase"]
path = "/Users/osbornejr/git/heterogeneous-graphlet-counting/packages/GraphletCounting"
uuid = "7ac45bc0-02f1-46da-ad35-65e91b15b4e1"
version = "0.1.0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "1cf1d7dcb4bc32d7b4a5add4232db3750c27ecb4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.8.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hexagons]]
deps = ["Test"]
git-tree-sha1 = "de4a6f9e7c4710ced6838ca906f81905f7385fd6"
uuid = "a1b4810d-1bce-5fbd-ac56-80944d57a21f"
version = "0.2.0"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "432b5b03176f8182bd6841fbfc42c718506a2d5f"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.15"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Infiltrator]]
deps = ["InteractiveUtils", "Markdown", "REPL", "UUIDs"]
git-tree-sha1 = "6e48065ac352c8c9616013faa419b0ea65bb6455"
uuid = "5903a43b-9cc3-4c30-8d17-598619ec4e9b"
version = "1.6.3"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.InformationMeasures]]
deps = ["Discretizers"]
git-tree-sha1 = "874d48f2026e8faf3fd55c86973fd028b02cd1a0"
uuid = "96684042-fbdc-5399-9b8e-d34e539a126c"
version = "0.3.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "Requires", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "42c17b18ced77ff0be65957a591d34f4ed57c631"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.31"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "6a125e6a4cb391e0b9adbd1afa9e771c2179f8ef"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.23"

[[deps.Juno]]
deps = ["Base64", "Logging", "Media", "Profile"]
git-tree-sha1 = "07cb43290a840908a771552911a6274bc6c072c7"
uuid = "e5e0dc1b-0480-54bc-9374-aad01c23163d"
version = "0.8.4"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "9816b296736292a80b9a3200eb7fbb57aaa3917a"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.5"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Librsvg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pango_jll", "Pkg", "gdk_pixbuf_jll"]
git-tree-sha1 = "ae0923dab7324e6bc980834f709c4cd83dd797ed"
uuid = "925c91fb-5dd6-59dd-8e8c-345e74382d89"
version = "2.54.5+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Loess]]
deps = ["Distances", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "46efcea75c890e5d820e670516dc156689851722"
uuid = "4345ca2d-374a-55d4-8d30-97f9976e7612"
version = "0.5.4"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "60168780555f3e663c536500aa790b6368adc02a"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.3.0"

[[deps.Luxor]]
deps = ["Base64", "Cairo", "Colors", "Dates", "FFMPEG", "FileIO", "Juno", "LaTeXStrings", "Random", "Requires", "Rsvg", "SnoopPrecompile"]
git-tree-sha1 = "909a67c53fddd216d5e986d804b26b1e3c82d66d"
uuid = "ae8d54c2-7ccd-5906-9d76-62fc9837b5bc"
version = "3.7.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "2ce8695e1e699b68702c03402672a69f54b8aca9"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Media]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "75a54abd10709c01f1b86b84ec225d26e840ed58"
uuid = "e89f7d12-3494-54d1-8411-f7d8b9ae1f27"
version = "0.5.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkConstruction]]
deps = ["Colors", "DataFrames", "DataPreprocessing", "Distributed", "Graphs", "InformationMeasures", "LinearAlgebra", "Luxor", "Printf", "ProgressMeter", "RCall", "SharedArrays", "Statistics", "StatsBase"]
path = "/Users/osbornejr/git/heterogeneous-graphlet-counting/packages/NetworkConstruction"
uuid = "6c2e41d2-72ae-425a-84e9-b8f08a301efb"
version = "0.1.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "82d7c9e310fe55aa54996e6f7f94674e2a38fcb4"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.9"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9ff31d101d987eb9d66bd8b176ac7c277beccd09"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.20+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "67eae2738d63117a196f497d7db789821bce61d1"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.17"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "84a314e3926ba9ec66ac097e3635e270986b0f10"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.9+0"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "5bb5129fdd62a2bbbe17c2756932259acf467386"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.50"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "LaTeXStrings", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "548793c7859e28ef026dba514752275ee871169f"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.ProjectFunctions]]
deps = ["CSV", "CategoricalArrays", "Colors", "Compose", "DataFrames", "DataPreprocessing", "DataStructures", "Dates", "Distributed", "Gadfly", "GraphPlot", "GraphletAnalysis", "GraphletCounting", "Graphs", "JLD2", "NetworkConstruction", "Pkg", "ProgressMeter", "Random", "StatsBase", "YAML"]
path = "/Users/osbornejr/git/heterogeneous-graphlet-counting/packages/ProjectFunctions"
uuid = "a8586eae-54f0-4952-9436-ba92c8ab3181"
version = "0.1.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "6ec7ac8412e83d57e313393220879ede1740f9ee"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.2"

[[deps.RCall]]
deps = ["CategoricalArrays", "Conda", "DataFrames", "DataStructures", "Dates", "Libdl", "Missings", "REPL", "Random", "Requires", "StatsModels", "WinReg"]
git-tree-sha1 = "d441bdeea943f8e8f293e0e3a78fe2d7c3aa24e6"
uuid = "6f49c342-dc21-5d91-9882-a32aef131414"
version = "0.13.15"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "feafdc70b2e6684314e188d95fe66d116de834a7"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.2"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.Rsvg]]
deps = ["Cairo", "Glib_jll", "Librsvg_jll"]
git-tree-sha1 = "3d3dc66eb46568fb3a5259034bfc752a0eb0c686"
uuid = "c4c386cf-5103-5370-be45-f3a111cca3b8"
version = "1.0.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "77d3c4726515dca71f6d80fbb5e251088defe305"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.18"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.ShiftedArrays]]
git-tree-sha1 = "503688b59397b3307443af35cd953a13e8005c16"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "2.0.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "63e84b7fdf5021026d0f17f76af7c57772313d99"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.21"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"

[[deps.StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "REPL", "ShiftedArrays", "SparseArrays", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "51cdf1afd9d78552e7a08536930d7abc3b288a5c"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.7.1"

[[deps.StringEncodings]]
deps = ["Libiconv_jll"]
git-tree-sha1 = "33c0da881af3248dafefb939a21694b97cfece76"
uuid = "69024149-9ee7-55f6-a4c4-859efe599b68"
version = "0.3.6"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "0b829474fed270a4b0ab07117dce9b9a2fa7581a"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.12"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WinReg]]
git-tree-sha1 = "cd910906b099402bcc50b3eafa9634244e5ec83b"
uuid = "1b915085-20d7-51cf-bf83-8f477d6f5128"
version = "1.0.0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.YAML]]
deps = ["Base64", "Dates", "Printf", "StringEncodings"]
git-tree-sha1 = "dbc7f1c0012a69486af79c8bcdb31be820670ba2"
uuid = "ddb6d928-2868-570f-bddf-ab3f9cf99eb6"
version = "0.4.8"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.gdk_pixbuf_jll]]
deps = ["Artifacts", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg", "Xorg_libX11_jll", "libpng_jll"]
git-tree-sha1 = "e9190f9fb03f9c3b15b9fb0c380b0d57a3c8ea39"
uuid = "da03df04-f53b-5353-a52f-6a8b0620ced0"
version = "2.42.8+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"
"""

# ╔═╡ Cell order:
# ╠═74c92bd0-ef4f-41a0-bc4b-5063ffcd12df
# ╠═6ea01822-2aa0-4073-84b0-304e5a86f9ea
# ╠═940181d4-a9b0-47e4-a13d-db2eb175e22e
# ╠═312f20e3-3c97-4108-847d-f2276abfab9a
# ╠═7211b91e-3ca1-4688-b88f-e9e5e92ca918
# ╠═72b2183e-99a9-4c9b-b74b-7e7966eb4bb8
# ╠═28836056-6604-4918-9f74-39bf81ad0559
# ╠═67f24598-4b27-4dd8-b533-45c1cc4e44a6
# ╠═4861b8cf-cf59-4cfd-b412-089ccfc00e90
# ╠═df21e8f4-7af6-4401-a8d7-e4f21e9c7139
# ╠═4ef32b7a-d7f3-4642-a75c-2870a5130ade
# ╠═8267a9ad-9551-4842-9194-5e250e79305e
# ╠═2733aa29-9e63-4b60-a7ee-0e9abcd62974
# ╠═e01833de-63e0-493a-8d6c-0240c1fe1633
# ╠═3c79ec7e-7e64-4f24-aa65-35c79167301b
# ╠═08edac1b-f869-4464-828d-bcb60c64a98e
# ╠═752f5e2a-6a63-442f-9e28-90db65b6eeb1
# ╠═6a2d3d61-cc8e-48ba-8839-705979033484
# ╠═812abebd-a67d-498c-939b-98f9b4d8a300
# ╠═b0ce6943-a214-4f5c-a75a-c3d59c71aa81
# ╠═0c721876-736e-495e-85b6-490ea7ded8f2
# ╠═2abc718c-e711-4b7f-880c-38922e225dd2
# ╠═4ee8ef83-ee2a-458c-b7ce-bf16fd4c5baa
# ╠═f7b11d20-9e7e-4197-b9b4-136d3e86644f
# ╠═822164ee-0193-4c16-a245-b55ea8082530
# ╠═18dbd2ba-ad4c-40fc-92bc-3b8579d8d952
# ╠═2c0de2b0-3d95-4abd-a7fe-7bbcf3fcfd52
# ╠═51ee4885-582a-4349-8d91-98f583b8a510
# ╠═8cdaf845-43c6-4455-822f-d30fbd2ff19b
# ╠═07904aca-d5da-444a-8fa5-55f6af49eb23
# ╠═7a0255c6-82c6-48c6-83e1-cd83606b05dc
# ╠═7f35a6f0-71e3-4448-8dd1-fe359ac73f94
# ╠═84aa8941-8766-4482-ba97-cfb5c477a072
# ╠═ebaa8810-9096-4cc2-9c28-31143143bf14
# ╠═3f3e6d45-d576-4385-bea8-e55a37d34512
# ╠═f34235b1-0e47-4f28-94db-ce3cfa598a91
# ╠═9bb259ec-8528-4049-80bf-5fa0d543e47c
# ╠═c2ddfc1b-5ef7-4c0f-82fb-d86c015aba74
# ╠═f24286b1-a29f-49dc-8005-058dfcf4440f
# ╠═c5f26162-6c29-46c2-b08b-8832ea301207
# ╠═33367e44-3233-4200-a091-a4bd4ea9adeb
# ╠═92c7926c-f21f-4665-84df-1ea98229cd4d
# ╠═534b73a8-ced9-4e8c-a9df-ef73af97198d
# ╠═c2d2fe61-9e00-4db8-aa80-97deb551486a
# ╠═87c1ea4a-26b2-42af-a21d-b6ff80366562
# ╠═3d0da8ba-1b04-462a-82e3-8b50eb1c29d8
# ╠═61a6f729-141c-4892-b983-dd48d78b6b7d
# ╠═726ac732-b9c4-439a-867c-61fcc34abedb
# ╠═f520c106-6e0d-49ec-8b9d-d656daa5b8dd
# ╠═ea2da746-ac59-4a54-99c1-5ab39a2f3221
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
