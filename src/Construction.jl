using Distributed, JLD, LightGraphs, GraphPlot, Colors, Random, Glob

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
	cache_dir = "$cwd/output/cache/$(params.test_name)_$(params.expression_cutoff)_$(params.norm_method)_$(params.variance_percent)_$(params.coexpression)_$(params.threshold)_$(params.threshold_method)"
	run(`mkdir -p $(cache_dir)`)
	## Check if there are cached files for runs with similar parameters (i.e. the same before irelevant parameters are invoked for a specific output)
	#similarity matrix
	sim_check = glob("output/cache/$(params.test_name)_$(params.expression_cutoff)_$(params.norm_method)_$(params.variance_percent)_$(params.coexpression)*/similarity*",cwd)
 	#check if any already exist
	if(length(sim_check)>0)
		#check if actual cache already exists
		if (!(cache_dir in first.(splitdir.(sim_check))))
			run(`cp $(sim_check[1]) $(cache_dir)/`)
		end
	end
	#Processing data:
	@info "Processing raw data..."
	raw_data = Array(select(raw_counts,filter(x->occursin("data",x),names(raw_counts))))
	## Clean - remove transcripts with total counts across all samples less than Cut
	##plot before cut
	histogram(DataFrame([log2.(vec(sum(raw_data,dims=2))),raw_counts[!,:transcript_type]],[:sum,:transcript_type]),:sum,:transcript_type,"$(params.website_dir)/_assets/$(params.page_name)/raw_data_histogram.svg",xaxis =" sum of expression (log2 adjusted)")
	
	clean_counts=raw_counts[vec(sum(raw_data,dims = 2 ).>=params.expression_cutoff),:]
	clean_data = Array(select(clean_counts,filter(x->occursin("data",x),names(clean_counts))))
	
	##plot after cut
	histogram(DataFrame([log2.(vec(sum(clean_data,dims=2))),clean_counts[!,:transcript_type]],[:sum,:transcript_type]),:sum,:transcript_type,"$(params.website_dir)/_assets/$(params.page_name)/clean_data_cut_histogram.svg",xaxis =" sum of expression (log2 adjusted)")
	
	#boxplot(raw_counts,"raw_data_cleaned_boxplot.svg")
	
	### Normalisation
	norm_data=library_size_normalisation(clean_data,params.norm_method)
	norm_counts = copy(clean_counts)
	norm_counts[:,findall(x->occursin("data",x),names(norm_counts))] = norm_data
	
	##Sampling for most variable transcripts
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
	adj_file = "$cache_dir/adjacency_matrix.jld"
	if (isfile(adj_file))
		@info "Loading adjacency matrix from $cache_dir..."
		pre_adj_matrix = JLD.load(adj_file,"pre-adj_matrix")
		adj_matrix = JLD.load(adj_file,"adjacency_matrix")
	else
		@info "Generating adjacency matrix..."
		if (params.threshold_method=="empirical_dist")
	 		pre_adj_matrix = empirical_dist_adjacency(similarity_matrix,params.threshold)
		if (params.threshold_method=="hard")
	 		pre_adj_matrix = adjacency(similarity_matrix,params.threshold)
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
	 	draw(SVG("$(params.website_dir)/_assets/alt_hub_network.svg",16cm,16cm),gplot(g,nodefillc = nodefillc))
	end	
	## Community structure
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
	
	if(params.graphlet_counting==true)

		graphlet_file = "$cache_dir/graphlets.jld" 
		if (isfile(graphlet_file))
			@info "Loading graphlet counts from $cache_dir..."
			graphlet_counts = JLD.load(graphlet_file,"graphlets")
		else
			@info "Counting graphlets..."
			@time graphlet_counts,Chi = count_graphlets(vertexlist,edgelist,4,run_method="distributed")
			#graphlet_concentrations = concentrate(graphlet_counts) 
			@info "Saving graphlet counts at $cache_dir..."
			##save the per-edge array as well in case we need it in the future (exp for debugging)
			JLD.save(graphlet_file,"graphlets",graphlet_counts,"Chi",Chi)
	
		end
		
	
		@info "Looking at typed representations of graphlets..."
		## randomise node types
		
		#number of randomised graphs
		N=100
		rand_types_set = [copy(vertexlist) for i in 1:N]
		#randomise each graph by node
		broadcast(shuffle!,rand_types_set) 
		
		rand_graphlets_file = "$cache_dir/rand_graphlets_$N.jld"
		if (isfile(rand_graphlets_file))
			@info "Loading randomised graphlet counts from $cache_dir..."
			rand_graphlet_collection = JLD.load(rand_graphlets_file,"rand graphlets")
		else
			@info "Counting graphlets on null model" 
			rand_graphlet_counts = count_graphlets.(rand_types_set,Ref(edgelist),4,run_method="distributed")
			rand_graphlet_dicts = broadcast(first,rand_graphlet_counts)
			rand_graphlet_collection = vcat(collect.(rand_graphlet_dicts)...)
			JLD.save(rand_graphlets_file,"rand graphlets",rand_graphlet_collection)
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
		for (i,hog) in enumerate(hom_graphlets)
			hog_df= DataFrame()
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
				histogram(rand_fil_fil,:log_value,:graphlet,"$(params.website_dir)/_assets/$(params.page_name)/plots/$(heg)_$(hog)_histogram.svg")
				rand_vals = rand_fil_fil[!,:value]
				rand_exp = sum(rand_vals)/N
				real_obs = real_dict[heg*"_"*hog]
				append!(hog_df,DataFrame(Graphlet = heg*"_"*hog, Expected = rand_exp,Observed = real_obs))	
			end
			##take log values to plot
			log_real_fil = copy(real_fil)
			log_real_fil.value =log.(log_real_fil.value)
			log_rand_fil = copy(rand_fil)
			log_rand_fil.value =log.(log_rand_fil.value)
			p = plot(layer(filter(:graphlet=>x->occursin(hog,x),log_real_fil),x = :graphlet,y = :value, Geom.point,color=["count in graph"]),Guide.xticks(label=true),Theme(key_position = :none),Guide.xlabel(nothing),Guide.ylabel("log value"),Guide.yticks(orientation=:vertical),layer(filter(:graphlet=>x->occursin(hog,x),log_rand_fil),x=:graphlet,y=:value,Geom.boxplot(suppress_outliers = true),color=:graphlet));
			draw(SVG("$(params.website_dir)/_assets/$(params.page_name)/plots/$(hog)_boxplot.svg",4inch,6inch),p)
			hog_array[i] = hog_df
		end
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
		
		#@time motif_counts = find_motifs(edgelist,"triangle_edge",100, typed = true, typelist = vec(vertexlist),plotfile="test.svg",graphlet_size = 4)
	end
end
