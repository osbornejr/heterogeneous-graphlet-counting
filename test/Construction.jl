### Include all source files TODO make this occur more fluently and automatically by creating a package, and using Revise
for src in filter(x->endswith(x,".jl"),readdir("src"))
    include(ENV["JULIA_PROJECT"]*"/src/"*src)
end

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
#create data matrix
raw_data=Array(raw_counts[!,2:25]);

#boxplot(raw_counts,"raw_data_boxplot.svg")

## Condense - merge polyA- and polyA+ counts for the same sample
condensed=raw_counts[!,[1:13;26]]
for i in 2:13
    condensed[!,i]=raw_counts[!,i]+raw_counts[!,i+12]
end
raw_counts=condensed
raw_data=Array(raw_counts[!,2:13]);

### Clean and normalise data- generate (filtered) raw_counts, norm_counts and normalised data matrix 
## Clean - remove transcripts with total counts across all samples less than X
X = 25
clean_counts=raw_counts[vec(sum(raw_data,dims = 2 ).>=X),:]
clean_data=Array(clean_counts[!,2:13]);


### Normalisation
norm_method = "median"
norm_data=library_size_normalisation(clean_data,norm_method)

## update data to normalised version
#data=per_sample
norm_counts=copy(clean_counts)
norm_counts[!,2:13]=norm_data


### Sample data - generate sample set of transcripts norm_data_sample, alongside reduced data matrix
#add variance column to normalised data
variance = vec(var(norm_data, dims=2))
insertcols!(norm_counts,"variance"=>variance)

#percent to select
X = 0.01
sample_counts_noncoding=sort(norm_counts[norm_counts[:transcript_type].=="noncoding",:],:variance)[Int(round(end*(1-X))):end,:]
sample_counts_coding=sort(norm_counts[norm_counts[:transcript_type].=="coding",:],:variance)[Int(round(end*(1-X))):end,:]
sample_counts = outerjoin(sample_counts_noncoding,sample_counts_coding,on = names(norm_counts))
sample_data=Array(sample_counts[!,2:13])


##Measure of coexpression
#similarity_matrix=mutual_information(data)
## store similarity matrices for use later:
coexpression = "mutual_information"
similarity_matrix = coexpression_measure(sample_data,coexpression)
## Adjacency matrix (using empricial distribution method atm)
threshold = 0.95
threshold_method = "empirical_dist"
if (threshold_method=="empirical_dist")
    pre_adj_matrix = empirical_dist_adjacency(similarity_matrix,threshold)
end
if (threshold_method=="hard")
    pre_adj_matrix = adjacency(similarity_matrix,threshold)
end
## Adjacency matrix

#maintain list of vertices in graph
vertexlist = sample_counts[:transcript_type]
##Consensus measure method
#adj_matrix = consensus_measure(data,methods = ["pearson","spearman","kendall","mutual_information"])

#Trim nodes with degree zero
network_counts = sample_counts[vec(sum(pre_adj_matrix,dims=2).!=0),:]
network_data = sample_data[vec(sum(pre_adj_matrix,dims=2).!=0),:]

vertexlist = vertexlist[vec(sum(pre_adj_matrix,dims=2).!=0)]
adj_matrix = copy(pre_adj_matrix)
adj_matrix = adj_matrix[:,vec(sum(pre_adj_matrix,dims=1).!=0)]
adj_matrix = adj_matrix[vec(sum(pre_adj_matrix,dims=2).!=0),:]
edgelist = edgelist_from_adj(adj_matrix)


###Sanity check... run the motif detection on an ALREADY random set of edges on the same node set
#sane_size = length(edgelist) 
#vertex_type_list = vec(vertexlist[:,2])
#rand1 = rand(1:size(vertexlist,1),sane_size)
#rand2 = rand(1:size(vertexlist,1),sane_size)
##first make sure there are no self loops
#sane_edgelist = splat(Pair).(eachrow(hcat(rand1[BitArray(((rand1.==rand2).-1).*-1)],rand2[BitArray(((rand1.==rand2).-1).*-1)])))
#
#now we order edges so that lower vertex label is first
#for (i,p) in enumerate(sane_edgelist)
#   if (first(p)>last(p))
#           sane_edgelist[i] = Pair(last(p),first(p))
#           end
#end
##finally, get rid of multiedges
#sane_edgelist = unique(sane_edgelist)
#

#get largest connected component TODO maybe have a way to select all components above a certain size, and plot them all separately?
#using LightGraphs
#using PrettyTables
g = Graph(adj_matrix) 
components = connected_components(g)
largest = components[length.(components).==max(length.(components)...)]
adj_matrix_comp = adj_matrix[largest[1],largest[1]]
g_comp = Graph(adj_matrix_comp)
##update vertexlist
vertexlist_comp = vertexlist[largest[1]]
##Network visualisation
#edgelist_comp=edgelist_from_adj(adj_matrix_comp)
#cytoscape_elements(vertexlist_comp,edgelist_comp,"cytoscape/elements.js")
#
#Network visualisation
using LightGraphs, GraphPlot
#g = SimpleGraph(adj_matrix)

#network plot
nodefillc = [colorant"lightseagreen", colorant"orange"][(vertexlist_comp.=="coding").+1]
draw(SVG("Mayank_$(norm_method)_$(threshold_method)_$(X)_$(coexpression)_component.svg",16cm,16cm),gplot(g_comp,nodefillc = nodefillc))
nodefillc = [colorant"lightseagreen", colorant"orange"][(vertexlist.=="coding").+1]
draw(SVG("Mayank_$(norm_method)_$(threshold_method)_$(X)_$(coexpression).svg",16cm,16cm),gplot(g,nodefillc = nodefillc))
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
@time graphlet_counts = count_graphlets(vertexlist,edgelist,4,run_method="distributed")
#graphlet_concentrations = concentrate(graphlet_counts) 

@time motif_counts = find_motifs(sane_edgelist,"hetero_rewire",1, typed = true, typelist = vertexlist,plotfile="test.svg",graphlet_size = 4)

