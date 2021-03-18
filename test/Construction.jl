### Include all source files TODO make this occur more fluently and automatically by creating a package, and using Revise
for src in filter(x->endswith(x,".jl"),readdir("src"))
	include(ENV["JULIA_PROJECT"]*"/src/"*src)
end

### Read in each set- generate (condensed) raw_counts and data matrix 
include(ENV["JULIA_PROJECT"]*"/test/ReadData.jl")

### Clean and normalise data- generate (filtered) raw_counts, norm_counts and normalised data matrix 
include(ENV["JULIA_PROJECT"]*"/test/CleanData.jl")

### Sample data - generate sample set of transcripts norm_data_sample, alongside reduced data matrix
include(ENV["JULIA_PROJECT"]*"/test/SampleData.jl")


#maintain list of vertices in graph
vertexlist = Array(norm_counts_sample[:,[1,14]])
##Measure of coexpression
#similarity_matrix=mutual_information(data)
similarity_matrix = coexpression_measure(data,"pearson")
## Adjacency matrix
threshold = 0.95
adj_matrix = adjacency(similarity_matrix,threshold)

##Consensus measure method
#adj_matrix = consensus_measure(data,methods = ["pearson","spearman","kendall","mutual_information"])


#Trim nodes with degree zero
data = data[vec(sum(adj_matrix,dims=2).!=0),:]
vertexlist = vertexlist[vec(sum(adj_matrix,dims=2).!=0),:]
adj_matrix = adj_matrix[:,vec(sum(adj_matrix,dims=1).!=0)]
adj_matrix = adj_matrix[vec(sum(adj_matrix,dims=2).!=0),:]
edgelist = edgelist_from_adj(adj_matrix)


##Sanity check... run the motif detection on an ALREADY random set of edges on the same node set
sane_size = length(edgelist) 
vertex_type_list = vec(vertexlist[:,2])
rand1 = rand(1:size(vertexlist,1),sane_size)
rand2 = rand(1:size(vertexlist,1),sane_size)
#first make sure there are no self loops
sane_edgelist = splat(Pair).(eachrow(hcat(rand1[BitArray(((rand1.==rand2).-1).*-1)],rand2[BitArray(((rand1.==rand2).-1).*-1)])))

#now we order edges so that lower vertex label is first
for (i,p) in enumerate(sane_edgelist)
	if (first(p)>last(p))
           sane_edgelist[i] = Pair(last(p),first(p))
       	end
end
#finally, get rid of multiedges
sane_edgelist = unique(sane_edgelist)
#get largest connected component TODO maybe have a way to select all components above a certain size, and plot them all separately?
#using LightGraphs
#using PrettyTables
#g = Graph(adj_matrix) 
#components = connected_components(g)
#largest = components[length.(components).==max(length.(components)...)]
#adj_matrix_comp = adj_matrix[largest[1],largest[1]]
#g_comp = Graph(adj_matrix_comp)
##update vertexlist
#vertexlist_comp = vertexlist[largest[1],:]
##Network visualisation
#edgelist_comp=edgelist_from_adj(adj_matrix_comp)
#cytoscape_elements(vertexlist_comp,edgelist_comp,"cytoscape/elements.js")
#
#Network Analysis
degrees = sum(adj_matrix,dims=2)
#p = plot(DataFrame(sort(degrees,dims=1)),x = "x1",Geom.histogram,Guide.title("Degree Distribution"),Guide.xlabel("degree"));
#draw(SVG("degree_distribution.svg"),p)
#connected_components_html_table(adj_matrix,"cytoscape/connected_components.html")
#
#Graphlet counting

##set up for distributed mode
#first clean to make sure there are no stray workers already around
using Distributed
rmprocs(workers())
#add workers equal to the number of available cpus	
addprocs(Threads.nthreads())
#addprocs(8)
@everywhere include(ENV["JULIA_PROJECT"]*"/src/GraphletCounting.jl")
@everywhere include(ENV["JULIA_PROJECT"]*"/src/NullModel.jl")
@time graphlet_counts = count_graphlets(vertexlist[:,2],edgelist,4,"distributed")
#graphlet_concentrations = concentrate(graphlet_counts) 

@time motif_counts = find_motifs(sane_edgelist,"hetero_rewire",1, typed = true, typelist = vertexlist[:,2],plotfile="test.svg",graphlet_size = 4)

