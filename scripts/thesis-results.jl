using Pkg
Pkg.activate(".")
using ProjectFunctions
using Results
experiment = "mayank-merged-1400-network"
#experiment = "GSE68559_sub"
raw_counts,round_counts,vst_counts,clean_counts,norm_counts,processed_counts = get_preprocessed_data("config/run-files/$(experiment).yaml")
Results.variance_histogram(norm_counts)

components,adj_matrix,network_counts,vertexlist,edgelist = get_network_construction()
fig = Results.plot_network(adj_matrix,vertex_colors = replace(vertexlist,"noncoding"=>:blue,"coding"=>:purple))

f = Results.typed_degree_distribution(vertexlist,edgelist)

#Results.add_to_fig(fig)
