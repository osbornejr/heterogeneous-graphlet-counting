for src in filter(x->endswith(x,".jl"),readdir("src"))
include("src/"*src)
end

### Read in each set
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
data=Array(raw_counts[!,2:25]);

boxplot(raw_counts,"raw_data_boxplot.svg")

test = stack(raw_counts,variable_name = "sample")

## Condense - merge polyA- and polyA+ counts for the same sample
condensed=raw_counts[!,[1:13;26]]
for i in 2:13
	condensed[!,i]=raw_counts[!,i]+raw_counts[!,i+12]
end
raw_counts=condensed
data=Array(raw_counts[!,2:13]);

boxplot(raw_counts,"raw_data_condensed_boxplot.svg")

## Clean - remove transcripts with total counts across all samples less than X
X = 25
raw_counts=raw_counts[vec(sum(data,dims = 2 ).>=X),:]
data=Array(raw_counts[!,2:13]);

boxplot(raw_counts,"raw_data_cleaned_boxplot.svg")

### Normalisation
data=library_size_normalisation(data,"upperquartile")
PCs,D_1,per_sample=pca(data')
#p=pca_plot(PCs,3);
#draw(SVG("output/Construction/test_per_sample.svg"),p)

## update data to normalised version
#data=per_sample
norm_counts=copy(raw_counts)
norm_counts[!,2:13]=data

boxplot(norm_counts,"norm_data_boxplot.svg")

p = plot(x = log2.(vec(data)), Geom.histogram(bincount = 100,density = false),Guide.xlabel("normalised expression (log2)"),Guide.ylabel("frequency"),Guide.title("Normalised expression histogram"));
draw(SVG("norm_expression.svg"),p)

norm_means = sum(data,dims=2)./12

p = plot(x = log2.(vec(norm_means)), Geom.histogram(bincount = 100,density = false),Guide.xlabel("normalised mean expression (log2)"),Guide.ylabel("frequency"),Guide.title("Normalised mean expression histogram"));
draw(SVG("norm_mean_expression.svg"),p)

## There might be some problems with keeping the polyA+ ad polyA- samples separate to this point. Monitor.


#To run locally (ie on laptop) we must still cut down number of transcripts to managable level here. TODO put into function? At the moment, we select the N/m highest variable transcripts of each type, where m is the number of types. 
N=10000

variance = vec(var(data, dims=2))
insertcols!(norm_counts,"variance"=>variance)

p = plot(x = log2.(variance), Geom.histogram(bincount = 100,density = false),Guide.xlabel("variance (log2)"),Guide.ylabel("frequency"),Guide.title("Variance histogram"));
draw(SVG("variance.svg"),p)

norm_counts_sample_noncoding=sort(norm_counts[norm_counts[:transcript_type].=="noncoding",:],:variance)[end-Int(N/2)+1:end,:]
norm_counts_sample_coding=sort(norm_counts[norm_counts[:transcript_type].=="coding",:],:variance)[end-Int(N/2)+1:end,:]
norm_counts_sample = outerjoin(norm_counts_sample_noncoding,norm_counts_sample_coding,on = names(norm_counts))
data=Array(norm_counts_sample[!,2:13])
#maintain list of vertices in graph
vertexlist = Array(norm_counts_sample[:,[1,14]])

##Measure of coexpression
#similarity_matrix=mutual_information(data)
similarity_matrix = coexpression_measure(data,"mutual_information")
## Adjacency matrix
threshold = 0.95
adj_matrix = adjacency(similarity_matrix,threshold)

##Consensus measure method
#adj_matrix = consensus_measure(data,methods = ["pearson","spearman","kendall","mutual_information"])

#Trim nodes with degree zero
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

@time graphlet_counts = count_graphlets(vertexlist[:,2],edgelist)
graphlet_concentrations = concentrate(graphlet_counts) 

#@time find_motifs(adj_matrix,"hetero_rewire",100, typed = true, typelist = vertexlist[:,2],plotfile="test.svg")

