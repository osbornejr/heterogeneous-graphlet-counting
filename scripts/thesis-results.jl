using Pkg
Pkg.activate(".")
using ProjectFunctions
using GraphletCounting
using StatsBase
using Results
using GLMakie
##set experiment
experiment = "mayank-merged-1400-network"
#experiment = "mayank-unmerged"
#experiment = "mayank-merged-altered-network"
#experiment = "mayank-merged"
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
#largest component
l_comp = components[findall(==(maximum(length.(components))),length.(components))...]
id_cut = 100
high_ids = l_comp[dd.>id_cut]
ids = l_comp[dd.<id_cut]
#compare processed_data for high degree to rest of proto-network
pd = data_from_dataframe(processed_counts)

## at this stage lets rearrange for the visualisation to have the sample columns in order (baked in order is 10,11,12,1,2,3,4,5,6,7,8,9; this translates to ICCV2-Stressed,JG11-Control,JG11-Stressed,ICCV2-Control. we will move the first three columns to the back to get each strain next to its condition pair.
pd = pd[:,[4,5,6,7,8,9,10,11,12,1,2,3]]
sample_names = ["JG11-Control-R1","JG11-Control-R2","JG11-Control-R3","JG11-Stress-R1","JG11-Stress-R2","JG11-Stress-R3","ICCV2-Control-R1","ICCV2-Control-R2","ICCV2-Control-R3","ICCV2-Stress-R1","ICCV2-Stress-R2","ICCV2-Stress-R3"]

f2 = Figure()
ax1 = Axis(f2[1,1],xticks=((2:3:12),["Control","Stress","Control","Stress"]),title="Expression counts: high degree transcripts")
#add extra layer of xtick labels
Axis(f2[1,1],xticks=((2:6:12),["JG11 (salt tolerant)","ICCV2 (salt sensitive)"]),xticklabelpad=25)
heatmap!(ax1,pd[high_ids,:]')
ax2 = Axis(f2[1,2],xticks=((2:3:12),["Control","Stress","Control","Stress"]),title="Expression counts: other transcripts")
Axis(f2[1,2],xticks=((2:6:12),["JG11 (salt tolerant)","ICCV2 (salt sensitive)"]),xticklabelpad=25)
heatmap!(ax2,pd[ids,:]')
f2

# it is clear that PCIT is not working as it should, and it appears the issue is the parital correlation step. If we try to get partialcor matrix manually in julia if also fails, as the underlying cor matrix is not PSD and thus not invertible. The underlying issue seems to be too many highly correlated transcripts. Lets identify them here.

#take correlation matrix
cor_pd = cor(pd')
#find all transcript 9pairs with cor greater than threshold
cor_threshold = 0.99
matching_pairs = filter(t->t[1]!=t[2],Tuple.(findall(>(cor_threshold),cor_pd)))
high_cors = unique(vcat(first.(matching_pairs),last.(matching_pairs)))
ax3 = Axis(f2[2,1],xticks=((2:3:12),["Control","Stress","Control","Stress"]),title="Expression counts: highly correlated (>$(cor_threshold)) transcripts")
Axis(f2[2,1],xticks=((2:6:12),["JG11 (salt tolerant)","ICCV2 (salt sensitive)"]),xticklabelpad=25)
heatmap!(ax3,pd[high_cors,:]')

##observe those transcripts that have high cor interactions but not high degree 
non_high_deg_high_cors = high_cors[.!in.(high_cors,Ref(intersect(high_ids,high_cors)))]
ax4 = Axis(f2[2,2],xticks=((2:3:12),["Control","Stress","Control","Stress"]),title="Expression counts: highly correlated but not high degree transcripts")
Axis(f2[2,2],xticks=((2:6:12),["JG11 (salt tolerant)","ICCV2 (salt sensitive)"]),xticklabelpad=25)
heatmap!(ax4,pd[non_high_deg_high_cors,:]')

non_high_cor_high_deg = high_ids[map(x->!in(x,intersect(high_ids,high_cors)),high_ids)]
ax5 = Axis(f2[3,1],xticks=((2:3:12),["Control","Stress","Control","Stress"]),title="Expression counts: high degree but not high cor transcripts")
Axis(f2[3,1],xticks=((2:6:12),["JG11 (salt tolerant)","ICCV2 (salt sensitive)"]),xticklabelpad=25)
heatmap!(ax5,pd[non_high_cor_high_deg,:]')

##find those transcripts that have an expression pattern that is on or off depending on strain (with a tolerance of 1 sample not following pattern)
#matching_strain_pattern = collect(1:size(pd)[1])[abs.(sum((pd.==0.0)[:,1:6],dims=2)-sum((pd.==0.0)[:,7:12],dims=2)).>4]
matching_JG_strain_pattern = collect(1:size(pd)[1])[sum((pd.==0.0)[:,1:6],dims=2)-sum((pd.==0.0)[:,7:12],dims=2).>4]
matching_ICCV_strain_pattern = collect(1:size(pd)[1])[sum((pd.==0.0)[:,1:6],dims=2)-sum((pd.==0.0)[:,7:12],dims=2).<-4]
matching_strain_pattern = vcat(matching_JG_strain_pattern,matching_ICCV_strain_pattern)
ax6 = Axis(f2[3,2],xticks=((2:3:12),["Control","Stress","Control","Stress"]),title="Expression counts: expression profiles matching strain samples")
Axis(f2[3,2],xticks=((2:6:12),["JG11 (salt tolerant)","ICCV2 (salt sensitive)"]),xticklabelpad=25)
heatmap!(ax6,pd[matching_strain_pattern,:]')
