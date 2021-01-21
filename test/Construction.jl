### Include all source files TODO make this occur more fluently and automatically by creating a package, and using Revise
for src in filter(x->endswith(x,".jl"),readdir("src"))
include("src/"*src)
end

### Read in each set- generate (condensed) raw_counts and data matrix 
include("test/ReadData.jl")

### Clean and normalise data- generate (filtered) raw_counts, norm_counts and normalised data matrix 
include("test/CleanData.jl")

### Sample data - generate sample set of transcripts norm_data_sample, alongside reduced data matrix
include("test/SampleData.jl")


#maintain list of vertices in graph
vertexlist = Array(norm_counts_sample[:,[1,14]])
##Measure of coexpression
#similarity_matrix=mutual_information(data)
similarity_matrix = coexpression_measure(data,"pcit")
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

#get largest connected component TODO maybe have a way to select all components above a certain size, and plot them all separately?
using LightGraphs
#using PrettyTables
g = Graph(adj_matrix) 
components = connected_components(g)
largest = components[length.(components).==max(length.(components)...)]
adj_matrix_comp = adj_matrix[largest[1],largest[1]]
g_comp = Graph(adj_matrix_comp)
#update vertexlist
vertexlist_comp = vertexlist[largest[1],:]
#Network visualisation
edgelist_comp=edgelist_from_adj(adj_matrix_comp)
cytoscape_elements(vertexlist_comp,edgelist_comp,"cytoscape/elements.js")

#Network Analysis
degrees = sum(adj_matrix,dims=2)
p = plot(DataFrame(sort(degrees,dims=1)),x = "x1",Geom.histogram,Guide.title("Degree Distribution"),Guide.xlabel("degree"));
draw(SVG("degree_distribution.svg"),p)
connected_components_html_table(adj_matrix,"cytoscape/connected_components.html")

#Graphlet counting

@time graphlet_counts = count_graphlets(vertexlist[:,2],edgelist,4)
#graphlet_concentrations = concentrate(graphlet_counts) 

#@time motif_counts = find_motifs(adj_matrix,"hetero_rewire",100, typed = true, typelist = vertexlist[:,2],plotfile="test.svg")

