### Include all source files TODO make this occur more fluently and automatically by creating a package, and using Revise
cwd = ENV["JULIA_PROJECT"]
for src in filter(x->endswith(x,".jl"),readdir("src"))
	include("$cwd/src/"*src)
end
#set up output directory
test_name = "GSE68559"
run(`mkdir -p "$cwd/output/$test_name"`)
using JLD


##set up for distributed mode
#first clean to make sure there are no stray workers already around
using Distributed
rmprocs(workers())
#add workers equal to the number of available cpus	
addprocs(Threads.nthreads())
#addprocs(8)
@everywhere include("$cwd/src/CoexpressionMeasures.jl")
@everywhere include("$cwd/src/GraphletCounting.jl")
@everywhere include("$cwd/src/NullModel.jl")


#Read in raw counts (cached)
raw_counts_file = "$cwd/output/cache/$(test_name)_raw_counts.jld"
if (isfile(raw_counts_file))
	raw_counts = JLD.load(raw_counts_file,"raw counts")
else
 	samples = CSV.read.(filter(x->occursin(".txt",x),readdir("$cwd/data/GSE68559_RAW",join=true)))
	sample_names = replace.(filter(x->occursin(".txt",x),readdir("$cwd/data/GSE68559_RAW")),"_isoforms_expr.txt"=>"").*" data"
	transcript_names = select(samples[1],1)
	raw_counts = rename!(hcat(transcript_names,select.(samples,Symbol("FPKM"))...,makeunique=true),["transcript_id";sample_names])
	raw_counts.transcript_type = replace(x-> occursin("lnc",x) ? "noncoding" : "coding",raw_counts.transcript_id)
	JLD.save(raw_counts_file,"raw counts",raw_counts)
end
raw_data = Array(select(raw_counts,filter(x->occursin("data",x),names(raw_counts))))


##plot before cut
histogram(DataFrame([log2.(vec(sum(raw_data,dims=2))),raw_counts[!,:transcript_type]],[:sum,:transcript_type]),:sum,:transcript_type,"output/$(test_name)/raw_data_histogram.svg",xaxis =" sum of expression (log2 adjusted)")

## Clean - remove transcripts with total counts across all samples less than Cut
Cut = 25
clean_counts=raw_counts[vec(sum(raw_data,dims = 2 ).>=Cut),:]
clean_data = Array(select(clean_counts,filter(x->occursin("data",x),names(clean_counts))))

##plot after cut
histogram(DataFrame([log2.(vec(sum(clean_data,dims=2))),clean_counts[!,:transcript_type]],[:sum,:transcript_type]),:sum,:transcript_type,"output/$(test_name)/clean_data_$(Cut)_cut_histogram.svg",xaxis =" sum of expression (log2 adjusted)")

#boxplot(raw_counts,"raw_data_cleaned_boxplot.svg")

### Normalisation
norm_method = "upper_quartile"
norm_data=library_size_normalisation(clean_data,norm_method)
norm_counts = copy(clean_counts)
norm_counts[:,findall(x->occursin("data",x),names(norm_counts))] = norm_data

##Sampling for most variable transcripts
#add variance column to normalised data
variance = vec(var(norm_data, dims=2))
norm_counts.variance = variance
X = 0.01
sample_counts_noncoding=sort(norm_counts[norm_counts[:transcript_type].=="noncoding",:],:variance)[Int(round(end*(1-X))):end,:]
sample_counts_coding=sort(norm_counts[norm_counts[:transcript_type].=="coding",:],:variance)[Int(round(end*(1-X))):end,:]
sample_counts = outerjoin(sample_counts_noncoding,sample_counts_coding,on = names(norm_counts))
sample_data = Array(select(sample_counts,filter(x->occursin("data",x),names(sample_counts))))

##Network construction
##Measure of coexpression
#similarity_matrix=mutual_information(data)
## file to cache similarity matrix for use later:
coexpression = "PID"
sim_file = "$cwd/output/cache/$(test_name)_similarity_matrix_$(norm_method)_$(X)_$(coexpression).jld"
if (isfile(sim_file))
	similarity_matrix = JLD.load(sim_file,"$(coexpression)_similarity_matrix")
else
	similarity_matrix = coexpression_measure(sample_data,coexpression)
	JLD.save(sim_file,"$(coexpression)_similarity_matrix",similarity_matrix)
end

## Adjacency matrix (using empricial distribution method atm)
threshold = 0.95
threshold_method = "empirical_dist"
if (threshold_method=="empirical_dist")
 	pre_adj_matrix = empirical_dist_adjacency(similarity_matrix,threshold)
end
if (threshold_method=="hard")
 	pre_adj_matrix = adjacency(similarity_matrix,threshold)
end
#maintain list of vertices in graph
vertexlist = sample_counts[:transcript_type]
#Trim nodes with degree zero
network_counts = sample_counts[vec(sum(pre_adj_matrix,dims=2).!=0),:]
network_data = sample_data[vec(sum(pre_adj_matrix,dims=2).!=0),:]

vertexlist = vertexlist[vec(sum(pre_adj_matrix,dims=2).!=0)]
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
draw(SVG("$cwd/output/pages/_assets/$(test_name)_$(norm_method)_$(threshold_method)_$(X)_$(coexpression)_component_network.svg",16cm,16cm),gplot(g_comp,nodefillc = nodefillc))
nodefillc = [colorant"lightseagreen", colorant"orange"][(vertexlist.=="coding").+1]
draw(SVG("$cwd/output/pages/_assets/$(test_name)_$(norm_method)_$(threshold_method)_$(X)_$(coexpression)_network.svg",16cm,16cm),gplot(g,nodefillc = nodefillc))

#Network Analysis

#Degrees
#homogonous degree distribution
degrees = vec(sum(adj_matrix,dims=2))
p = plot(DataFrame([sort(degrees)]),x = "x1",Geom.histogram,Guide.title("Degree distribution"),Guide.xlabel("degree"));
draw(SVG("$cwd/output/pages/_assets/$(test_name)_degree_distribution.svg"),p)

#degrees for each transcript type
for type in unique(vertexlist)
	p = plot(DataFrame([sort(degrees[vertexlist.==type])]),x = "x1",Geom.histogram,Guide.title("Degree distribution"),Guide.xlabel("degree"));
	draw(SVG("$cwd/output/pages/_assets/$(test_name)_$(type)_degree_distribution.svg"),p)
end

## Hubs
deg_thresh = mean(degrees)+2*std(degrees)
nodefillc = [colorant"black", colorant"red"][(degrees.>deg_thresh).+1]
draw(SVG("$cwd/output/pages/_assets/$(test_name)_$(norm_method)_$(threshold_method)_$(X)_$(coexpression)_degree_network.svg",16cm,16cm),gplot(g,nodefillc = nodefillc))

##HTML output (concept atm)
run(`mkdir -p output/$(test_name)/page`)
open("output/$(test_name)/page/index.html","w") do io
	show(io,"text/html",gplot(g,nodefillc = nodefillc))
	show(io,"text/html",p)
end

@time graphlet_counts = count_graphlets(vertexlist,edgelist,4,run_method="distributed")
#graphlet_concentrations = concentrate(graphlet_counts) 


@time motif_counts = find_motifs(edgelist,"triangle_edge",100, typed = true, typelist = vec(vertexlist),plotfile="test.svg",graphlet_size = 4)
