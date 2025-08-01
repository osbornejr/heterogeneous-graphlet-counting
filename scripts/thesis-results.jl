using Pkg
Pkg.activate(".")
using ProjectFunctions
using Results
using GLMakie
##set experiment
#experiment = "mayank-merged-1400-network"
experiment = "mayank-merged"
#experiment = "GSE68559_sub"

##load preprocessing data
raw_counts,round_counts,vst_counts,clean_counts,norm_counts,processed_counts = get_preprocessed_data("config/run-files/$(experiment).yaml")
#plot variance histogram
Results.variance_histogram(norm_counts)


##load network info 
components,adj_matrix,network_counts,vertexlist,edgelist = get_network_construction()
#visualise network
fig = Results.plot_network(adj_matrix,vertex_colors = replace(vertexlist,"noncoding"=>:blue,"coding"=>:purple))
#view typed degree distributions
f = Results.typed_degree_distribution(vertexlist,edgelist)

#Results.add_to_fig(fig)
#


#finding high degree nodes
#degree distribution
dd = sum.(values.(GraphletCounting.typed_degree_distribution(vertexlist,edgelist)))
#get ids for large component
high_ids = components[1][dd.>100]
ids = components[1][dd.<100]
#compare processed_data for high degree to rest of proto-network
f2 = Figure()
ax1 = Axis(f2[1,1],title="Expression counts: high degree transcripts")
heatmap!(ax1,data_from_dataframe(processed_counts[high_ids,:])')
ax2 = Axis(f2[1,2],title="Expression counts: other transcripts")
heatmap!(ax2,data_from_dataframe(processed_counts[ids,:])')
f2
