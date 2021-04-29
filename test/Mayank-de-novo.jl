### Include all source files TODO make this occur more fluently and automatically by creating a package, and using Revise
cwd = ENV["JULIA_PROJECT"]
for src in filter(x->endswith(x,".jl"),readdir("src"))
	include("$cwd/src/"*src)
end
##Run parameters
test_name = "Mayank-de-novo"
site_name = "TestWebsite"
page_name = "menu2"
Cut = 25
norm_method = "upper_quartile"
variance_percent = 0.01
coexpression = "pcit"
threshold = 0.95
threshold_method = "empirical_dist"

run_parameter_df = DataFrame(Test = test_name, Expression_cut_off = Cut,Normalisation = norm_method, Variance_cut_off = variance_percent, Coexpression_measure = coexpression, Edge_threshold = threshold, Threshold_method = threshold_method)
CSV.write("$cwd/output/$(site_name)/_assets/$(page_name)/tableinput/run_parameters.csv",run_parameter_df)

#set up output directories
run(`mkdir -p "$cwd/output/$test_name"`)
run(`mkdir -p "$cwd/output/$(site_name)/_assets/$(page_name)/tableinput"`)
run(`mkdir -p "$cwd/output/$(site_name)/_assets/$(page_name)/plots"`)
using JLD


##set up for distributed mode
#first clean to make sure there are no stray workers already around
using Distributed
if(length(workers())!=Threads.nthreads())
	rmprocs(workers())
 	#add workers equal to the number of available cpus	
	addprocs(Threads.nthreads())
 	#addprocs(8)
	@everywhere include("src/CoexpressionMeasures.jl")
	@everywhere include("src/GraphletCounting.jl")
	@everywhere include("src/NullModel.jl")
end

#Read in raw counts (cached)
raw_counts_file = "$cwd/output/cache/$(test_name)_raw_counts.jld"
if (isfile(raw_counts_file))
	raw_counts = JLD.load(raw_counts_file,"raw counts")
else
	### Read in each set- generate (condensed) raw_counts and data matrix 
	raw_counts=read_count_data("data/mayank-de-novo/isoforms",method="expected_count");
	#raw_counts=RSEM.read_count_data("data/mayank-per-transcript/isoforms",method="expected_count");
	
	#filtering out into types
	code_counts=filter_count_data("data/mayank-de-novo/code-hits.list",raw_counts)
	noncode_counts=filter_count_data("data/mayank-de-novo/non-code-hits.list",raw_counts)	
	insertcols!(code_counts,"transcript_type"=>"coding")
	insertcols!(noncode_counts,"transcript_type"=>"noncoding")
	
	#merging back into one dataframe, with 
	raw_counts=outerjoin(code_counts,noncode_counts,on = intersect(names(code_counts),names(noncode_counts)))
	
	
	#boxplot(raw_counts,"raw_data_boxplot.svg")
	
	## Condense - merge polyA- and polyA+ counts for the same sample
	condensed=raw_counts[!,[1:13;26]]
	for i in 2:13
		condensed[!,i]=raw_counts[!,i]+raw_counts[!,i+12]
	end
	raw_counts=condensed
	JLD.save(raw_counts_file,"raw counts",raw_counts)
end
raw_data = Array(select(raw_counts,filter(x->occursin("data",x),names(raw_counts))))


##plot before cut
histogram(DataFrame([log2.(vec(sum(raw_data,dims=2))),raw_counts[!,:transcript_type]],[:sum,:transcript_type]),:sum,:transcript_type,"output/$(site_name)/_assets/$(page_name)/raw_data_histogram.svg",xaxis =" sum of expression (log2 adjusted)")

## Clean - remove transcripts with total counts across all samples less than Cut
clean_counts=raw_counts[vec(sum(raw_data,dims = 2 ).>=Cut),:]
clean_data = Array(select(clean_counts,filter(x->occursin("data",x),names(clean_counts))))

##plot after cut
histogram(DataFrame([log2.(vec(sum(clean_data,dims=2))),clean_counts[!,:transcript_type]],[:sum,:transcript_type]),:sum,:transcript_type,"output/$(site_name)/_assets/$(page_name)/clean_data_$(Cut)_cut_histogram.svg",xaxis =" sum of expression (log2 adjusted)")

#boxplot(raw_counts,"raw_data_cleaned_boxplot.svg")

### Normalisation
norm_data=library_size_normalisation(clean_data,norm_method)
norm_counts = copy(clean_counts)
norm_counts[:,findall(x->occursin("data",x),names(norm_counts))] = norm_data

##Sampling for most variable transcripts
#add variance column to normalised data
variance = vec(var(norm_data, dims=2))
norm_counts.variance = variance
sample_counts_noncoding=sort(norm_counts[norm_counts[:transcript_type].=="noncoding",:],:variance)[Int(round(end*(1-variance_percent))):end,:]
sample_counts_coding=sort(norm_counts[norm_counts[:transcript_type].=="coding",:],:variance)[Int(round(end*(1-variance_percent))):end,:]
sample_counts = outerjoin(sample_counts_noncoding,sample_counts_coding,on = names(norm_counts))
sample_data = Array(select(sample_counts,filter(x->occursin("data",x),names(sample_counts))))

##Network construction
##Measure of coexpression
#similarity_matrix=mutual_information(data)
## file to cache similarity matrix for use later:
sim_file = "$cwd/output/cache/$(test_name)_similarity_matrix_$(norm_method)_$(X)_$(coexpression).jld"
if (isfile(sim_file))
	similarity_matrix = JLD.load(sim_file,"$(coexpression)_similarity_matrix")
else
	similarity_matrix = coexpression_measure(sample_data,coexpression)
	JLD.save(sim_file,"$(coexpression)_similarity_matrix",similarity_matrix)
end

## Adjacency matrix (using empricial distribution method atm)
if (threshold_method=="empirical_dist")
 	pre_adj_matrix = empirical_dist_adjacency(similarity_matrix,threshold)
end
if (threshold_method=="hard")
 	pre_adj_matrix = adjacency(similarity_matrix,threshold)
end
#Trim nodes with degree zero
network_counts = sample_counts[vec(sum(pre_adj_matrix,dims=2).!=0),:]
network_data = sample_data[vec(sum(pre_adj_matrix,dims=2).!=0),:]
#maintain list of vertices in graph
vertex_names = network_counts[:transcript_id]
vertexlist = copy(network_counts[:transcript_type])

##form final adjacency matrix
adj_matrix = copy(pre_adj_matrix)
adj_matrix = adj_matrix[:,vec(sum(pre_adj_matrix,dims=1).!=0)]
adj_matrix = adj_matrix[vec(sum(pre_adj_matrix,dims=2).!=0),:]
edgelist = edgelist_from_adj(adj_matrix)

#Network visualisation
using LightGraphs, GraphPlot
g = SimpleGraph(adj_matrix)
##get largest component
components = connected_components(g)
largest = components[length.(components).==max(length.(components)...)]
adj_matrix_comp = adj_matrix[largest[1],largest[1]]
g_comp = Graph(adj_matrix_comp)
##update vertexlist
vertexlist_comp = vertexlist[largest[1]]
##plot (either connected component or whole network
nodefillc = [colorant"lightseagreen", colorant"orange"][(vertexlist_comp.=="coding").+1]
draw(SVG("$cwd/output/$(site_name)/_assets/$(page_name)/$(norm_method)_$(threshold_method)_$(X)_$(coexpression)_component_network.svg",16cm,16cm),gplot(g_comp,nodefillc = nodefillc))
nodefillc = [colorant"lightseagreen", colorant"orange"][(vertexlist.=="coding").+1]
draw(SVG("$cwd/output/$(site_name)/_assets/$(page_name)/$(norm_method)_$(threshold_method)_$(X)_$(coexpression)_network.svg",16cm,16cm),gplot(g,nodefillc = nodefillc))

#Network Analysis
#Type representations 
##set up csv string
csv = "Step,Coding counts,Non-coding counts,Non-coding proportion\n"
csv = csv*"Raw counts,"*string(size(raw_counts,1)-size(filter(:transcript_type=>x->x == "noncoding",raw_counts),1))*","*string(size(filter(:transcript_type=>x->x == "noncoding",raw_counts),1))*","*string(round(size(filter(:transcript_type=>x->x == "noncoding",raw_counts),1)/size(raw_counts,1),sigdigits=3))*"\n"
csv = csv*"Clean counts,"*string(size(clean_counts,1)-size(filter(:transcript_type=>x->x == "noncoding",clean_counts),1))*","*string(size(filter(:transcript_type=>x->x == "noncoding",clean_counts),1))*","*string(round(size(filter(:transcript_type=>x->x == "noncoding",clean_counts),1)/size(clean_counts,1),sigdigits=3))*"\n"
csv = csv*"Sample counts,"*string(size(sample_counts,1)-size(filter(:transcript_type=>x->x == "noncoding",sample_counts),1))*","*string(size(filter(:transcript_type=>x->x == "noncoding",sample_counts),1))*","*string(round(size(filter(:transcript_type=>x->x == "noncoding",sample_counts),1)/size(sample_counts,1),sigdigits=3))*"\n"
csv = csv*"Network counts,"*string(size(network_counts,1)-size(filter(:transcript_type=>x->x == "noncoding",network_counts),1))*","*string(size(filter(:transcript_type=>x->x == "noncoding",network_counts),1))*","*string(round(size(filter(:transcript_type=>x->x == "noncoding",network_counts),1)/size(network_counts,1),sigdigits=3))*"\n"
write("$cwd/output/$(site_name)/_assets/$(page_name)/tableinput/type_representation.csv",csv)

#Degrees
#homogonous degree distribution
degrees = vec(sum(adj_matrix,dims=2))
p = plot(DataFrame([sort(degrees)]),x = "x1",Geom.histogram,Guide.title("Degree distribution"),Guide.xlabel("degree"));
draw(SVG("$cwd/output/$(site_name)/_assets/$(page_name)/degree_distribution.svg"),p)

#degrees for each transcript type
for type in unique(vertexlist)
	p = plot(DataFrame([sort(degrees[vertexlist.==type])]),x = "x1",Geom.histogram,Guide.title("Degree distribution"),Guide.xlabel("degree"));
	draw(SVG("$cwd/output/$(site_name)/_assets/$(page_name)/$(type)_degree_distribution.svg"),p)
end

## Hubs
deg_thresh = mean(degrees)+2*std(degrees)
nodefillc = [colorant"black", colorant"red"][(degrees.>deg_thresh).+1]
draw(SVG("$cwd/output/$(site_name)/_assets/$(page_name)/$(norm_method)_$(threshold_method)_$(X)_$(coexpression)_degree_network.svg",16cm,16cm),gplot(g,nodefillc = nodefillc))
deg_thresh = 70#mean(degrees)+2*std(degrees)
nodefillc = [colorant"black", colorant"red"][(degrees.>deg_thresh).+1]
draw(SVG("$cwd/output/$(site_name)/_assets/$(page_name)/$(norm_method)_$(threshold_method)_$(X)_$(coexpression)_degree70_network.svg",16cm,16cm),gplot(g,nodefillc = nodefillc))

## Community structure
using RCall
@rput adj_matrix
@rput vertex_names
@rput site_name
@rput page_name
R"""
library(RColorBrewer)
library(igraph)
library(threejs)
#keep this last to mask other packages:
library(htmlwidgets)
library(biomaRt)
library(GO.db)
library(httr)
library(topGO)
library(tidyverse)

##change to tibble
vertex_names = tibble(vertex_names)
vertex_names_trimmed = sapply(vertex_names,tools::file_path_sans_ext)
pre_g <- graph.adjacency(adj_matrix,mode="undirected",diag=F)
#Now we add attributes to graph by rebuilding it via edge and vertex lists
edges <- as_tibble(as.data.frame(get.edgelist(pre_g)))
vertices <- tibble(name = 1:nrow(vertex_names),vertex_names)%>% dplyr::rename(label=vertex_names)

g <- graph_from_data_frame(edges,directed = F, vertices)
#
##delete zero degree vertices
g <- delete.vertices(g,igraph::degree(g)==0)
##find maximum connected component
g <- decompose(g, mode = "strong", max.comps = NA, min.vertices = 10)[[1]]
#
###Community structure
##get only connected vertices
vertices <- vertices %>% filter(name %in% names(V(g)))

##cluster graph and colour by community
communities <- cluster_louvain(g)

##add colours to graph
vertices <- vertices %>% mutate(group = as_factor(communities$membership),color = group)
##if there are more than 11 communities, spectral colour palette is not sufficient. so we concat two palettes (if there are more than 22 comms, need another!)
x = length(unique(communities$membership))-11
colour_palette = c(brewer.pal(name = "Spectral", n = 11),brewer.pal(name = "BrBG", n = x))
levels(vertices$color) <- colour_palette
edges <- as_tibble(as.data.frame(get.edgelist(g)))
g <- graph_from_data_frame(edges,directed = F, vertices)
vertex_attr(g,"size") <- 0.5
plot = graphjs(g,bg = "white");
saveWidget(plot,paste0("output/",site_name,"/",page_name,"_communities.html"))
##also save separate from website
saveWidget(plot,paste0("output/cache/",page_name,"_communities.html"))
## Functional annotation of communities
#get list of string arrays for transcript ids in each community
comms <- split(vertices,vertices$group)
comm_ids = lapply(comms,function(x) x$label)
comm_ids_trimmed = sapply(comm_ids,function(y) sapply(y,function(x) tools::file_path_sans_ext(x)))
"""
@rget comms
@rget vertices




@time graphlet_counts = count_graphlets(vertexlist,edgelist,4,run_method="distributed")
#graphlet_concentrations = concentrate(graphlet_counts) 

## randomise node types
using Random
#number of randomised graphs
N=100
rand_types_set = [copy(vertexlist) for i in 1:N]
#randomise each graph by node
broadcast(shuffle!,rand_types_set) 

rand_graphlets_file = "$cwd/output/cache/$(test_name)_rand_graphlets_$N.jld"
if (isfile(rand_graphlets_file))
	rand_graphlet_collection = JLD.load(rand_graphlets_file,"rand graphlets")
else
	rand_graphlet_counts = count_graphlets.(rand_types_set,Ref(edgelist),4,run_method="distributed")
	rand_graphlet_dicts = broadcast(first,rand_graphlet_counts)
	rand_graphlet_collection = vcat(collect.(rand_graphlet_dicts)...)
	JLD.save(rand_graphlets_file,"rand graphlets",rand_graphlet_collection)
end


rand_df = DataFrame(graphlet = broadcast(first,rand_graphlet_collection),value = broadcast(last,rand_graphlet_collection))
real_df = DataFrame(graphlet = broadcast(first,collect(graphlet_counts[1])),value = broadcast(last,collect(graphlet_counts[1])))

##function to get the n permutations of a set xs
all_perm(xs, n) = vec(map(collect, Iterators.product(ntuple(_ -> xs, n)...)))

##convert graphlet_counts dict output to default dictionary, returning 0 for graphlets that don't exist in the real network
real_dict = DefaultDict(0,graphlet_counts[1])
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
		
		rand_vals = filter(:graphlet=>x->x==heg*"_"*hog,rand_df)[!,:value]
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
	draw(SVG("$cwd/output/TestWebsite/_assets/$(page_name)/plots/$(hog)_histogram.svg",4inch,6inch),p)
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

@time motif_counts = find_motifs(edgelist,"triangle_edge",100, typed = true, typelist = vec(vertexlist),plotfile="test.svg",graphlet_size = 4)