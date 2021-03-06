using Distributed, JLD, LightGraphs, GraphPlot, Colors, Random, Glob, Distributions

function distributed_setup(inclusions::Array{String,1})

	##set up for distributed mode

	#first clean to make sure there are no stray workers already around
	if(length(workers())!=Threads.nthreads())
		rmprocs(workers())
	 	#add workers equal to the number of available cpus	
		addprocs(Threads.nthreads())
	 	#addprocs(8)
		@everywhere inclusions
		for inc in inclusions
			@everywhere include(inc)
		end
	end
end
function webpage_construction(raw_counts::DataFrame,params::RunParameters)


	@info "Building directory structure..."
	##establish output directories	
	run(`mkdir -p "$(params.website_dir)/_assets/$(params.page_name)/tableinput"`)
	run(`mkdir -p "$(params.website_dir)/_assets/$(params.page_name)/plots"`)

	run_parameter_df = DataFrame(Test = params.test_name, Expression_cut_off = params.expression_cutoff,Normalisation = params.norm_method, Variance_cut_off = params.variance_percent, Coexpression_measure = params.coexpression, Edge_threshold = params.threshold, Threshold_method = params.threshold_method)
	CSV.write("$(params.website_dir)/_assets/$(params.page_name)/tableinput/run_parameters.csv",run_parameter_df)
	#set up output directories
	run(`mkdir -p "$(params.website_dir)/$(params.page_name)"`)
	run(`mkdir -p "$(params.website_dir)/_assets/$(params.page_name)/tableinput"`)
	run(`mkdir -p "$(params.website_dir)/_assets/$(params.page_name)/plots"`)
		
	##Cache setup
	#New method: cache directory updates folder by folder as we go. Avoids need for moving files around,symbolic links etc. 
	cache_dir = "$cwd/output/cache/$(params.test_name)"
	run(`mkdir -p $(cache_dir)`)
	## Check if there are cached files for runs with similar parameters (i.e. the same before irelevant parameters are invoked for a specific output)
	#TODO improve this process to automatically detect via a dependency tree (will also tidy up cache dir)
	#similarity matrix
	#sim_check = glob("output/cache/$(params.test_name)_$(params.expression_cutoff)_$(params.norm_method)_$(params.variance_percent)_$(params.coexpression)*/similarity*",cwd)
 	#check if any already exist
	#if(length(sim_check)>0)
		#check if actual cache already exists
	#	if (!(cache_dir in first.(splitdir.(sim_check))))
	#		run(`cp $(sim_check[1]) $(cache_dir)/`)
	#	end
	#end
	#functional annotation
	#func_check = glob("output/cache/$(params.test_name)_$(params.expression_cutoff)_$(params.norm_method)_$(params.variance_percent)_$(params.coexpression)_$(params.threshold)_$(params.threshold_method)*/simil_$(params.threshold)_$(params.threshold_method)*/func*",cwd)
 	#check if any already exist
	#if(length(func_check)>0)
		#check if actual cache already exists
	#	if (!(cache_dir in first.(splitdir.(func_check))))
	#		run(`cp $(func_check[1]) $(cache_dir)/`)
	#	end
	#end
	#Processing data:
	@info "Processing raw data..."
	raw_data = Array(select(raw_counts,filter(x->occursin("data",x),names(raw_counts))))
	
	## Clean - remove transcripts with total counts across all samples less than Cut
	##update cache
	cache_dir = "$cwd/output/cache/$(params.test_name)/cutoff/$(params.expression_cutoff)"
	run(`mkdir -p $(cache_dir)`)
	##plot before cut
	histogram(DataFrame([log2.(vec(sum(raw_data,dims=2))),raw_counts[!,:transcript_type]],[:sum,:transcript_type]),:sum,:transcript_type,"$(params.website_dir)/_assets/$(params.page_name)/raw_data_histogram.svg",xaxis =" sum of expression (log2 adjusted)")
	
	clean_counts=raw_counts[vec(sum(raw_data,dims = 2 ).>=params.expression_cutoff),:]
	clean_data = Array(select(clean_counts,filter(x->occursin("data",x),names(clean_counts))))
	
	##plot after cut
	histogram(DataFrame([log2.(vec(sum(clean_data,dims=2))),clean_counts[!,:transcript_type]],[:sum,:transcript_type]),:sum,:transcript_type,"$(params.website_dir)/_assets/$(params.page_name)/clean_data_cut_histogram.svg",xaxis =" sum of expression (log2 adjusted)")
	
	#boxplot(raw_counts,"raw_data_cleaned_boxplot.svg")
	
	### Normalisation
	##update cache
	cache_dir = "$cwd/output/cache/$(params.test_name)/cutoff/$(params.expression_cutoff)/normalisation/$(params.norm_method)"
	run(`mkdir -p $(cache_dir)`)
	norm_data=library_size_normalisation(clean_data,params.norm_method)
	norm_counts = copy(clean_counts)
	norm_counts[:,findall(x->occursin("data",x),names(norm_counts))] = norm_data
	
	##Sampling for most variable transcripts
	##update cache
	cache_dir = "$cwd/output/cache/$(params.test_name)/cutoff/$(params.expression_cutoff)/normalisation/$(params.norm_method)/sampling/$(params.variance_percent)"
	run(`mkdir -p $(cache_dir)`)
	#add variance column to normalised data
	variance = vec(var(norm_data, dims=2))
	norm_counts.variance = variance
	sample_counts_noncoding=sort(norm_counts[norm_counts[:transcript_type].=="noncoding",:],:variance)[Int(round(end*(1-params.variance_percent))):end,:]
	sample_counts_coding=sort(norm_counts[norm_counts[:transcript_type].=="coding",:],:variance)[Int(round(end*(1-params.variance_percent))):end,:]
	sample_counts = outerjoin(sample_counts_noncoding,sample_counts_coding,on = names(norm_counts))
	sample_data = Array(select(sample_counts,filter(x->occursin("data",x),names(sample_counts))))
	##Network construction
	@info "Constructing the network..."
	##Measure of coexpression
	##update cache
	cache_dir = "$cwd/output/cache/$(params.test_name)/cutoff/$(params.expression_cutoff)/normalisation/$(params.norm_method)/sampling/$(params.variance_percent)/similarity/$(params.coexpression)"
	run(`mkdir -p $(cache_dir)`)
	#similarity_matrix=mutual_information(data)
	## file to cache similarity matrix for use later:
	sim_file = "$cache_dir/similarity_matrix.jld"
	if (isfile(sim_file))
		@info "Loading similarity matrix from $cache_dir..."
		similarity_matrix = JLD.load(sim_file,"similarity_matrix")
	else
		@info "Generating similarity matrix..."
		similarity_matrix = coexpression_measure(sample_data,params.coexpression)
		@info "Saving similarity matrix at $cache_dir..."
		JLD.save(sim_file,"similarity_matrix",similarity_matrix)
	end
	
	## Adjacency matrix 
	##update cache
	cache_dir = "$cwd/output/cache/$(params.test_name)/cutoff/$(params.expression_cutoff)/normalisation/$(params.norm_method)/sampling/$(params.variance_percent)/similarity/$(params.coexpression)/threshold/$(params.threshold)/threshold_method/$(params.threshold_method)"
	run(`mkdir -p $(cache_dir)`)
	adj_file = "$cache_dir/adjacency_matrix.jld"
	if (isfile(adj_file))
		@info "Loading adjacency matrix from $cache_dir..."
		pre_adj_matrix = JLD.load(adj_file,"pre-adj_matrix")
		adj_matrix = JLD.load(adj_file,"adjacency_matrix")
	else
		@info "Generating adjacency matrix..."
		if (params.threshold_method=="empirical_dist")
	 		pre_adj_matrix = empirical_dist_adjacency(similarity_matrix,params.threshold)
		elseif (params.threshold_method=="empirical_dist_zero")
	 		pre_adj_matrix = empirical_dist_zero_adjacency(similarity_matrix,params.threshold)
		elseif (params.threshold_method=="hard")
	 		pre_adj_matrix = adjacency(similarity_matrix,params.threshold)
		elseif (params.threshold_method=="top")
			##TODO setting top x value here for now; should be a parameter, but as an Int rather than Float as threshold param is for other methods
			pre_adj_matrix = top_adjacency(similarity_matrix,10)
		end
		##form final adjacency matrix
		adj_matrix = copy(pre_adj_matrix)
		adj_matrix = adj_matrix[:,vec(sum(pre_adj_matrix,dims=1).!=0)]
		adj_matrix = adj_matrix[vec(sum(pre_adj_matrix,dims=2).!=0),:]
		@info "Saving adjacency matrix at $cache_dir..."
		JLD.save(adj_file,"adjacency_matrix",adj_matrix,"pre-adj_matrix",pre_adj_matrix)
	end
	#Trim nodes with degree zero
	network_counts = sample_counts[vec(sum(pre_adj_matrix,dims=2).!=0),:]
	network_data = sample_data[vec(sum(pre_adj_matrix,dims=2).!=0),:]
	#maintain list of vertices in graph
	vertex_names = network_counts[:transcript_id]
	vertexlist = copy(network_counts[:transcript_type])	
	edgelist = edgelist_from_adj(adj_matrix)
	
	#Network visualisation
	g = SimpleGraph(adj_matrix)
	##get largest component
	components = connected_components(g)
	largest = components[length.(components).==max(length.(components)...)]
	adj_matrix_comp = adj_matrix[largest[1],largest[1]]
	g_comp = Graph(adj_matrix_comp)
	##update vertexlist
	vertexlist_comp = vertexlist[largest[1]]

	if (params.visualise==true)
		@info "Visualising network..."
		##plot (either connected component or whole network
		nodefillc = [colorant"lightseagreen", colorant"orange"][(vertexlist_comp.=="coding").+1]
		draw(SVG("$(params.website_dir)/_assets/$(params.page_name)/component_network.svg",16cm,16cm),gplot(g_comp,nodefillc = nodefillc))
		nodefillc = [colorant"lightseagreen", colorant"orange"][(vertexlist.=="coding").+1]
		draw(SVG("$(params.website_dir)/_assets/$(params.page_name)/network.svg",16cm,16cm),gplot(g,nodefillc = nodefillc))
	end
	#Network Analysis
	@info "Analysing network..."
	#Type representations 
	##set up csv string
	csv = "Step,Coding counts,Non-coding counts,Non-coding proportion\n"
	csv = csv*"Raw counts,"*string(size(raw_counts,1)-size(filter(:transcript_type=>x->x=="noncoding",raw_counts),1))*","*string(size(filter(:transcript_type=>x->x=="noncoding",raw_counts),1))*","*string(round(size(filter(:transcript_type=>x->x=="noncoding",raw_counts),1)/size(raw_counts,1),sigdigits=3))*"\n"
	csv = csv*"Clean counts,"*string(size(clean_counts,1)-size(filter(:transcript_type=>x->x=="noncoding",clean_counts),1))*","*string(size(filter(:transcript_type=>x->x=="noncoding",clean_counts),1))*","*string(round(size(filter(:transcript_type=>x->x=="noncoding",clean_counts),1)/size(clean_counts,1),sigdigits=3))*"\n"
	csv = csv*"Sample counts,"*string(size(sample_counts,1)-size(filter(:transcript_type=>x->x=="noncoding",sample_counts),1))*","*string(size(filter(:transcript_type=>x->x=="noncoding",sample_counts),1))*","*string(round(size(filter(:transcript_type=>x->x=="noncoding",sample_counts),1)/size(sample_counts,1),sigdigits=3))*"\n"
	csv = csv*"Network counts,"*string(size(network_counts,1)-size(filter(:transcript_type=>x->x=="noncoding",network_counts),1))*","*string(size(filter(:transcript_type=>x->x=="noncoding",network_counts),1))*","*string(round(size(filter(:transcript_type=>x->x=="noncoding",network_counts),1)/size(network_counts,1),sigdigits=3))*"\n"
	write("$(params.website_dir)/_assets/$(params.page_name)/tableinput/type_representation.csv",csv)
	
	#Degrees
	#homogonous degree distribution
	degrees = vec(sum(adj_matrix,dims=2))
	p = plot(DataFrame([sort(degrees)]),x = "x1",Geom.histogram,Guide.title("Degree distribution"),Guide.xlabel("degree"));
	draw(SVG("$(params.website_dir)/_assets/$(params.page_name)/degree_distribution.svg"),p)
	
	#degrees for each transcript type
	for type in unique(vertexlist)
		p = plot(DataFrame([sort(degrees[vertexlist.==type])]),x = "x1",Geom.histogram,Guide.title("Degree distribution"),Guide.xlabel("degree"));
		draw(SVG("$(params.website_dir)/_assets/$(params.page_name)/$(type)_degree_distribution.svg"),p)
	end
	
	## Hubs
	if (params.visualise==true)
		@info "Identifying hubs..."
	 	deg_thresh = Int(floor(mean(degrees)+2*std(degrees)))
	 	nodefillc = [colorant"black", colorant"red"][(degrees.>deg_thresh).+1]
	 	draw(SVG("$(params.website_dir)/_assets/$(params.page_name)/two_std_hub_network.svg",16cm,16cm),gplot(g,nodefillc = nodefillc))
	 	deg_thresh = 70#mean(degrees)+2*std(degrees)
	 	nodefillc = [colorant"black", colorant"red"][(degrees.>deg_thresh).+1]
		draw(SVG("$(params.website_dir)/_assets/$(params.page_name)/alt_hub_network.svg",16cm,16cm),gplot(g,nodefillc = nodefillc))
	end	
	## Community structure
	##update cache
	cache_dir = "$cwd/output/cache/$(params.test_name)/cutoff/$(params.expression_cutoff)/normalisation/$(params.norm_method)/sampling/$(params.variance_percent)/similarity/$(params.coexpression)/threshold/$(params.threshold)/threshold_method/$(params.threshold_method)/communities/"
	run(`mkdir -p $(cache_dir)`)
	@info "Identifying communities..."
	##use gene ids here, as they have more chance of getting a GO annotation
	if(params.func_annotate==true)
	 	vertex_gene_names = network_counts[:gene_id]
		community_vertices = get_community_structure(adj_matrix,vertex_gene_names,"louvain",threejs_plot = true,plot_prefix = "$(params.website_dir)/$(params.page_name)") 
		## functional annotations of communities
		func_file = "$cache_dir/func_annotations.jld" 
		if (isfile(func_file))
			functional_annotations = JLD.load(func_file,"functional annotations")
		else
			functional_annotations = get_functional_annotations(community_vertices,ensembl_version = "75",write_csv = true, csv_dir ="$(params.website_dir)/_assets/$(params.page_name)/tableinput/")	
			JLD.save(func_file,"functional annotations",functional_annotations)
		end
	else
		community_vertices = get_community_structure(adj_matrix,vertex_names,"louvain",threejs_plot = true,plot_prefix = "$(params.website_dir)/$(params.page_name)") 
	end
	
	##Network stats table
	##set up csv string
	csv = "Nodes,Edges,Components,Nodes in largest component,Maximal degree,communities detected\n"
	csv = csv*string(length(vertexlist))*","*string(length(edgelist))*","*string(length(components))*","*string(length(largest...))*","*string(max(degrees...))*","*string(length(unique(community_vertices.group)))*"\n"
	write("$(params.website_dir)/_assets/$(params.page_name)/tableinput/network_stats.csv",csv)
	if(params.graphlet_counting==true)

		##update cache
		cache_dir = "$cwd/output/cache/$(params.test_name)/cutoff/$(params.expression_cutoff)/normalisation/$(params.norm_method)/sampling/$(params.variance_percent)/similarity/$(params.coexpression)/threshold/$(params.threshold)/threshold_method/$(params.threshold_method)/communities/graphlets/"
		run(`mkdir -p $(cache_dir)`)
		graphlet_file = "$cache_dir/graphlets.jld" 
		if (isfile(graphlet_file))
			@info "Loading graphlet counts from $cache_dir..."
			graphlet_counts = JLD.load(graphlet_file,"graphlets")
			timer = JLD.load(graphlet_file,"time")
		else
			@info "Counting graphlets..."
			timer=@elapsed graphlet_counts = count_graphlets(vertexlist,edgelist,4,run_method="distributed-old")
			#graphlet_concentrations = concentrate(graphlet_counts) 
			@info "Saving graphlet counts at $cache_dir..."
			##save the per-edge array as well in case we need it in the future (exp for debugging)
			JLD.save(graphlet_file,"graphlets",graphlet_counts,"time",timer)
	
		end
	

		#method to deduce which run method the null model should use given the run time of graphlets above ($timer)
		if (timer<30)
			null_run = "distributed-short"
		else
			null_run = "distributed-long"
		end

 	 	##Typed representations
		@info "Looking at typed representations of graphlets..."
		
		##update cache
		cache_dir = "$cwd/output/cache/$(params.test_name)/cutoff/$(params.expression_cutoff)/normalisation/$(params.norm_method)/sampling/$(params.variance_percent)/similarity/$(params.coexpression)/threshold/$(params.threshold)/threshold_method/$(params.threshold_method)/communities/graphlets/typed_representations/nullmodel/$(params.null_model_size)_simulations"
		run(`mkdir -p $(cache_dir)`)
		N=params.null_model_size
		rand_graphlets_file = "$cache_dir/rand_graphlets.jld"
		if (isfile(rand_graphlets_file))
			@info "Loading randomised vertices and graphlet counts from $cache_dir..."
			rand_types_set = JLD.load(rand_graphlets_file,"rand vertices")
			rand_graphlet_collection = JLD.load(rand_graphlets_file,"rand graphlets")
		else
			## randomise node types
			
			#number of randomised graphs
			rand_types_set = [copy(vertexlist) for i in 1:N]
			#randomise each graph by node
			broadcast(shuffle!,rand_types_set) 
			@info "Counting graphlets on null model" 
			if (null_run=="distributed-short")
				#rand_graphlet_counts = count_graphlets.(rand_types_set,Ref(edgelist),4,run_method="distributed-old")
				rand_graphlet_counts = @showprogress pmap(x->count_graphlets(x,edgelist,4,run_method="serial"),rand_types_set,batch_size =10)
			end
			if (null_run=="distributed-long")
				#rand_graphlet_counts = count_graphlets.(rand_types_set,Ref(edgelist),4,run_method="distributed")
				rand_graphlet_counts = @showprogress map(x->count_graphlets(x,edgelist,4,run_method="distributed"),rand_types_set)
			end
			rand_graphlet_collection = vcat(collect.(rand_graphlet_counts)...)
			@info "Saving random graphlet count information at $cache_dir..."
			JLD.save(rand_graphlets_file,"rand graphlets",rand_graphlet_collection,"rand vertices",rand_types_set)
		end
		
		
		rand_df = DataFrame(graphlet = broadcast(first,rand_graphlet_collection),value = broadcast(last,rand_graphlet_collection))
		real_df = DataFrame(graphlet = broadcast(first,collect(graphlet_counts)),value = broadcast(last,collect(graphlet_counts)))
		
		##function to get the n permutations of a set xs
		all_perm(xs, n) = vec(map(collect, Iterators.product(ntuple(_ -> xs, n)...)))
		
		##convert graphlet_counts dict output to default dictionary, returning 0 for graphlets that don't exist in the real network
		real_dict = DefaultDict(0,graphlet_counts)
		hom_graphlets = unique(last.(split.(unique(real_df[:graphlet]),"_")))
		##array to store all homogonous graphlet dfs
		hog_array=Array{DataFrame,1}(undef,length(hom_graphlets))	
		hog_array_under=Array{DataFrame,1}(undef,length(hom_graphlets))	
		for (i,hog) in enumerate(hom_graphlets)
			hog_df= DataFrame()
			hog_df_under= DataFrame()
			## restrict info to just hg 
			real_fil = filter(:graphlet=>x->occursin(hog,x),real_df)
			rand_fil = filter(:graphlet=>x->occursin(hog,x),rand_df)
		
			#get hetero subgraphlets within homogonous type (problem: might not be complete set present in real/rand outputs?)
			if (occursin("4",hog))
				het_graphlets = union(first.(split.(real_fil[:graphlet],"_4")),first.(split.(rand_fil[:graphlet],"_4")))
			elseif (occursin("3",hog))
				het_graphlets = union(first.(split.(real_fil[:graphlet],"_3")),first.(split.(real_fil[:graphlet],"_3")))
			end 
			for heg in het_graphlets
				rand_fil_fil = filter(:graphlet=>x->x==heg*"_"*hog,rand_df)
				transform!(rand_fil_fil,:value =>ByRow(x-> log(x))=>:log_value)
				##histogram for each heterogeneous graphlet
				histogram(rand_fil_fil,:value,:graphlet,"$(params.website_dir)/_assets/$(params.page_name)/plots/$(heg)_$(hog)_histogram.svg")
				##log version
				histogram(rand_fil_fil,:log_value,:graphlet,"$(params.website_dir)/_assets/$(params.page_name)/plots/$(heg)_$(hog)_log_histogram.svg")
				rand_vals = rand_fil_fil[!,:value]
				log_rand_vals =rand_fil_fil[!,:log_value]
				rand_exp = sum(rand_vals)/N
				log_rand_exp = sum(log_rand_vals)/N
				real_obs = real_dict[heg*"_"*hog]
				log_real_obs = log(real_dict[heg*"_"*hog]+1)
				## Monte Carlo method
				r = sum(rand_vals.>=real_obs)
				r_under = sum(rand_vals.<=real_obs)
				##special case: real_obs is zero.. Zero occurences do not register in rand_vals table, but clearly every rand network has equal or greater count 
				if(real_obs==0)
					r = N
				end
				p_value = (r+1)/(N+1)
				p_value_under = (r_under+1)/(N+1)			
				##Z-score method (assumes either normal or lognormal distribution of counts in sims)
				#using real values
			#	z_score = (abs(real_obs) - rand_exp)/std(rand_vals)
				##using log values
				#z_score = (abs(log_real_obs) - log_rand_exp)/std(log_rand_vals)
				#p_value =pdf(Normal(0,1),z_score) 
				append!(hog_df,DataFrame(Graphlet = heg*"_"*hog, Expected = rand_exp,Observed = real_obs, p_value = p_value))	
				append!(hog_df_under,DataFrame(Graphlet = heg*"_"*hog, Expected = rand_exp,Observed = real_obs, p_value = p_value_under))	
			end
			##take log values to plot
			log_real_fil = copy(real_fil)
			log_real_fil.value =log.(log_real_fil.value)
			log_rand_fil = copy(rand_fil)
			log_rand_fil.value =log.(log_rand_fil.value)
			p = plot(layer(filter(:graphlet=>x->occursin(hog,x),log_real_fil),x = :graphlet,y = :value, Geom.point,color=["count in graph"]),Guide.xticks(label=true),Theme(key_position = :none),Guide.xlabel(nothing),Guide.ylabel("log value"),Guide.yticks(orientation=:vertical),layer(filter(:graphlet=>x->occursin(hog,x),log_rand_fil),x=:graphlet,y=:value,Geom.boxplot(suppress_outliers = true),color=:graphlet));
			draw(SVG("$(params.website_dir)/_assets/$(params.page_name)/plots/$(hog)_boxplot.svg",4inch,6inch),p)
			hog_array[i] = hog_df
			hog_array_under[i] = hog_df_under
		end
		## find significant graphlets
		sig_graphlets = vcat(filter.(:p_value=>x->x<0.05,hog_array)...)
		insig_graphlets = vcat(filter.(:p_value=>x->x<0.05,hog_array_under)...)
		
		#save in output cache
		html_table_maker(sig_graphlets,"$cache_dir/sig_type_representations.html",imgs=sig_graphlets.Graphlet)				
		html_table_maker(insig_graphlets,"$cache_dir/insig_type_representations.html",imgs=sig_graphlets.Graphlet)				
		#save for website version
		html_table_maker(sig_graphlets,"$(params.website_dir)/_assets/$(params.page_name)/sig_type_representations.html",imgs=sig_graphlets.Graphlet,figpath="../figs/")
		html_table_maker(insig_graphlets,"$(params.website_dir)/_assets/$(params.page_name)/insig_type_representations.html",imgs=insig_graphlets.Graphlet,figpath="../figs/")	
		##look at edge types in randomised networks
		real_type_edgecounts = countmap(splat(tuple).(sort.(eachrow(hcat(map(x->vertexlist[x],first.(edgelist)),map(x->vertexlist[x],last.(edgelist)))))))
		rand_types_edgecounts = map(y->(countmap(splat(tuple).(sort.(eachrow(hcat(map(x->y[x],first.(edgelist)),map(x->y[x],last.(edgelist)))))))),rand_types_set)
		rand_edge_collection = vcat(collect.(rand_types_edgecounts)...)
		rand_edge_df = DataFrame(graphlet = broadcast(first,rand_edge_collection),value = broadcast(last,rand_edge_collection))
		random_edges = DataFrame()
		for t in unique(rand_edge_df[:graphlet])
			rand_vals = filter(:graphlet=>x->x==t,rand_edge_df)[!,:value]
			rand_exp = sum(rand_vals)/N
			real_obs = real_type_edgecounts[t]
			append!(random_edges,DataFrame(Graphlet = first(t)*"_"*last(t)*"_edge", Expected = rand_exp,Observed = real_obs))	
		end
		
		#pretty_table(random_edges,backend=:html,standalone = false)
		
		#@time motif_counts = find_motifs(edgelist,"hetero_rewire",100, typed = true, typelist = vec(vertexlist),plotfile="$cache_dir/motif_detection.svg",graphlet_size = 4)

#		### Validation steps
#		# looking at identified significant graphlets and seeing if they check out biologically
#		#To do this we require the node-level relationships output of the graphlet counting algorithm. We it is generally quicker to run this than save and load it via JLD, as it is quite large) 
#		
#	#	if (isfile(graphlet_file))
#	#		@info "Loading graphlet counts from $cache_dir..."
#	#		Rel = JLD.load(graphlet_file,"Rel")
#	#	else
#	#		@info "Counting graphlets..."
#	#		##Do full relationships run in case per-node relationships are required (TODO add toggle for this?)
#	#		timer=@elapsed graphlet_counts,Chi,Rel = count_graphlets(vertexlist,edgelist,4,run_method="distributed-old",relationships = true)
#	#		#graphlet_concentrations = concentrate(graphlet_counts) 
#	#		@info "Saving graphlet counts at $cache_dir..."
#	#		##save the per-edge array as well in case we need it in the future (exp for debugging)
#	#		JLD.save(graphlet_file,"graphlets",graphlet_counts,"Chi",Chi,"Rel",Rel,"time",timer)
#	#
#	#	end
#
		#Now moved into a function
 	 	#Coincidents = get_KEGG_graphlet_coincidences(vertexlist,adj_matrix)
#
#
# 		graphlet_counts,Chi,Rel = count_graphlets(vertexlist,edgelist,4,run_method="distributed-old",relationships = true,progress = true)
# 		## combine relationships into one array
# 		rel = vcat(Rel...)
# 		#rel_array = broadcast(a->[i for i in a],broadcast(x->x[1:end-1],rel))
#		##remove 0 from 3-node entries
#		#rel_array = map(y->filter(x->x!=0,y),rel_array)
# 		#rel_names = map(x->broadcast(y->vertex_gene_names[y],x),rel_array)	
# 		#rel_transcript_names = map(x->broadcast(y->vertex_names[y],x),rel_array)	
#		
#
#		##split into each homogeneous graphlet and sort there
#		##sort and find unique copies of graphlet (TODO unique for each graphlet!)
# 		graphlet_types = string.(unique(last.(split.(collect(keys(graphlet_counts)),"_"))))
#		graphlet_rels = Dict{String,Array{Array{Int64,1},1}}()
#		@time for g in graphlet_types
# 			hogs = filter(x->x[end]==g,rel)
# 			hogs_array = broadcast(a->[i for i in a],broadcast(x->x[1:end-1],hogs))
#			if(g == "3-path")
#				#get rid of leading zero
#				hogs_array = map(y->filter(x->x!=0,y),hogs_array)
#				graphlet_rels[g] = unique(broadcast(x->[sort(x[[1,3]])[1],x[2],sort(x[[1,3]])[2]],hogs_array))
# 			end
#			if(g == "3-tri")
#				#get rid of leading zero
#				hogs_array = map(y->filter(x->x!=0,y),hogs_array)
#				graphlet_rels[g] = unique(broadcast(x->[sort(x[[1,2,3]])...],hogs_array))
# 			end
#			if(g == "4-path")
#				graphlet_rels[g] = unique(broadcast(x->[sort(x[[1,4]])[1],sort(x[[2,3]])...,sort(x[[1,4]])[2]],hogs_array))
# 			end
#			if(g == "4-star")
#				graphlet_rels[g] = unique(broadcast(x->[sort(x[[1,2,4]])[1:2]...,x[3],sort(x[[1,2,4]])[3]],hogs_array))
# 			end
#			if(g == "4-tail")
#				graphlet_rels[g] = unique(broadcast(x->[sort(x[[1,2]])...,x[3],x[4]],hogs_array))
# 			end
#			if(g == "4-cycle")
#				graphlet_rels[g] = unique(broadcast(x->[sort(x[[1,2,3,4]])...],hogs_array))
# 			end
#			if(g == "4-chord")
#				graphlet_rels[g] = unique(broadcast(x->[sort(x[[1,4]])[1],sort(x[[2,3]])...,sort(x[[1,4]])[2]],hogs_array))
# 			end
#			if(g == "4-clique")
#				graphlet_rels[g] = unique(broadcast(x->[sort(x[[1,2,3,4]])...],hogs_array))
# 			end
#		end
#		## add edges to graphlet_rels as separate entry
#		graphlet_rels["2-path"] = [[x...] for x in eachrow(hcat(first.(edgelist),last.(edgelist)))]
#
#		#graphlet of interest... set manually for now 
# 		goi = sig_graphlets.Graphlet[2]		
# 		hegoi,hogoi = string.(split(goi,"_")[1:end-1]),string(split(goi,"_")[end])
# 		#filter down to homogenous graphlet first
# 		hogs = filter(x->x[end]==hogoi,rel)
# 		hogs_array = broadcast(a->[i for i in a],broadcast(x->x[1:end-1],hogs))
# 		##sort and find unique copies of graphlet (TODO unique for each graphlet!)
# 		if(hogoi == "4-star")
# 			hogs_sorted = unique(broadcast(x->[sort(x[[1,2,4]])[1:2]...,x[3],sort(x[[1,2,4]])[3]],hogs_array))
# 		end
# 		
# 		#get those that match heterogeneous pattern as well
# 		hogs_types = map(x->broadcast(y->vertexlist[y],x),hogs_sorted)	
# 		hegs = hogs_sorted[findall(x->x== hegoi,hogs_types)]
# 	 	#get names of transcripts in matching pattern 
# 		hegs_names = map(x->broadcast(y->vertex_gene_names[y],x),hegs)	
# 		hegs_transcript_names = map(x->broadcast(y->vertex_names[y],x),hegs)	
#		#restart R session	
#		R"""
#		sapply(names(sessionInfo()$otherPkgs),function(pkg) detach(paste0('package:',pkg),character.only =T,force = T));rm(list=ls())
#		"""
#		#use small sample for now
#		hegs_sample = sample(hegs_names,20)
#		@rput hegs_sample
#		@rput vertex_names
#		@rput vertex_gene_names
#		R"""
#		library(biomaRt)
#		library(httr)
#		library(edgeR)
#		library(tidyverse)
#		hegs_sample_trimmed = lapply(hegs_sample,function (x) sapply(x,tools::file_path_sans_ext))
#		## connect to biomart
#		set_config(config(ssl_verifypeer = 0L))
#		ensembl_version = "current"	
#		if (ensembl_version=="current")
#			{
#			##mirrors to try: "useast" "uswest" "asia"
#			ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",mirror="uswest", dataset = "hsapiens_gene_ensembl") 
#			} else
#			{
#			ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",version=ensembl_version) 
#			}
#		hegs_ids <- lapply(hegs_sample_trimmed,function(x) getBM(attributes=c("ensembl_gene_id","entrezgene_id"),"ensembl_gene_id",x, mart = ensembl,useCache=FALSE)) 
#		keggs = lapply(hegs_ids,function(x) topKEGG(kegga(x[[2]]),n=15,truncate =34))
#		
#
#		
#
#		##Alternative approach: select KEGG terms of interest from whole set of transcripts first, and then look for those in graphlets
#		## full list of entrez_ids mapped to transcripts/genes (neither will be one-to-one, use whichever offers better coverage. Transcripts also preferred as a better reflection of the data rather than expanding to the gene level).
#		##transcripts
#		transcripts_trimmed = sapply(vertex_names,tools::file_path_sans_ext)	
#		transcripts = data.frame(trimmed_names = transcripts_trimmed)
#		entrez_from_transcripts = getBM(attributes=c("ensembl_transcript_id","ensembl_gene_id","entrezgene_id"),"ensembl_transcript_id",transcripts_trimmed, mart = ensembl,useCache=FALSE)	
#		transcript_coverage = length(entrez_from_transcripts[[3]])-sum(is.na(entrez_from_transcripts[[3]]))
#		transcripts$entrez_id = entrez_from_transcripts[match(transcripts_trimmed,entrez_from_transcripts[[1]]),3]
#		
#		##genes
#		genes_trimmed = sapply(vertex_gene_names,tools::file_path_sans_ext)	
#		genes = data.frame(trimmed_names = genes_trimmed)
#		entrez_from_genes = getBM(attributes=c("ensembl_gene_id","entrezgene_id"),"ensembl_gene_id",genes_trimmed, mart = ensembl,useCache=FALSE)	
#		gene_coverage = length(entrez_from_genes[[2]])-sum(is.na(entrez_from_genes[[2]]))
#		genes$entrez_id = entrez_from_genes[match(genes_trimmed,entrez_from_genes[[1]]),2]
#		
#
#		#get list of entrez ids mapped to KEGG pathways 
#		KEGG <-getGeneKEGGLinks(species.KEGG="hsa")
#		## get top hits to select from
#		top_terms = topKEGG(kegga(entrez_from_transcripts[[3]],n=Inf,truncate = 34))
#		##get the network candidates for each pathway
#		per_pathway = sapply(1:nrow(top_terms),function(x) KEGG$GeneID[ KEGG$PathwayID == row.names(top_terms)[x]])
#		in_network = lapply(test,function(x) transcripts$entrez_id %in% x)
#			
#
#		"""
#		@rget in_network
#		candidates = Dict{String,Array{Int,1}}()
#		for e in in_network
#			candidates[string(first(e))] = findall(.==(true),last(e))
#		end
#	 	
#		Coincidents = Dict{String,Dict{String,Array{Tuple,1}}}()
#		for ent in candidates
#			@info "checking candidates for $(first(ent))..."
#			cands = last(ent)
#			#one_coincidents = Array{Array{Int64,1}}(undef,length(keys(graphlet_rels))) 
#			#two_coincidents = Array{Array{Int64,1}}(undef,length(keys(graphlet_rels))) 
#			#three_coincidents = Array{Array{Int64,1}}(undef,length(keys(graphlet_rels))) 
#			#four_coincidents = Array{Array{Int64,1}}(undef,length(keys(graphlet_rels))) 
#			#one_coincidents = Array{Tuple,1}()
#			two_coincidents = Array{Tuple,1}()
#			three_coincidents = Array{Tuple,1}()
#			four_coincidents = Array{Tuple,1}()
#			for g in keys(graphlet_rels) 
#				#graphlets with at least two candidate transcripts involved
#				#one_coincidents[i] = findall(x->sum(map(y->in(y,x),cands))>0,graphlet_rels[g])
#				#two_coincidents[i] = findall(x->sum(map(y->in(y,x),cands))>1,graphlet_rels[g])
#				#three_coincidents[i] = findall(x->sum(map(y->in(y,x),cands))>2,graphlet_rels[g])
#				#four_coincidents[i] = findall(x->sum(map(y->in(y,x),cands))>3,graphlet_rels[g])
#				#push!(one_coincidents,map(x->tuple(graphlet_rels[g][x]...,g),findall(x->sum(map(y->in(y,x),cands))>0,graphlet_rels[g])...))
#				push!(two_coincidents,map(x->tuple(graphlet_rels[g][x]...,g),findall(x->sum(map(y->in(y,x),cands))>1,graphlet_rels[g]))...)
#				push!(three_coincidents,map(x->tuple(graphlet_rels[g][x]...,g),findall(x->sum(map(y->in(y,x),cands))>2,graphlet_rels[g]))...)
#				push!(four_coincidents,map(x->tuple(graphlet_rels[g][x]...,g),findall(x->sum(map(y->in(y,x),cands))>3,graphlet_rels[g]))...)
#			end
#			#save each in dictionary for candidate
#			Coincidents[first(ent)] = Dict("two"=>two_coincidents,"three"=>three_coincidents,"four"=>four_coincidents) 
#		end
		#@time motif_counts = find_motifs(edgelist,"hetero_rewire",100, typed = true, typelist = vec(vertexlist),plotfile="$cache_dir/motif_detection.svg",graphlet_size = 4)
	end
end
