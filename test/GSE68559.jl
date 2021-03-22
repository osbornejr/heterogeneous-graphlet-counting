### Include all source files TODO make this occur more fluently and automatically by creating a package, and using Revise
for src in filter(x->endswith(x,".jl"),readdir("src"))
	include(ENV["JULIA_PROJECT"]*"/src/"*src)
end
samples = CSV.read.(filter(x->occursin(".txt",x),readdir("data/GSE68559_RAW",join=true)))
sample_names = replace.(filter(x->occursin(".txt",x),readdir("data/GSE68559_RAW")),"_isoforms_expr.txt"=>"").*" data"
transcript_names=select(samples[1],1)
raw_counts=rename!(hcat(transcript_names,select.(samples,Symbol("FPKM"))...,makeunique=true),["transcript_id";sample_names])

raw_counts.transcript_type = replace(x-> occursin("lnc",x) ? "noncoding" : "coding",raw_counts.transcript_id)
raw_data = Array(select(raw_counts,filter(x->occursin("data",x),names(raw_counts))))


## Clean - remove transcripts with total counts across all samples less than X
X = 25
clean_counts=raw_counts[vec(sum(raw_data,dims = 2 ).>=X),:]
clean_data = Array(select(norm_counts,filter(x->occursin("data",x),names(norm_counts))))

#boxplot(raw_counts,"raw_data_cleaned_boxplot.svg")

### Normalisation
norm_data=library_size_normalisation(clean_data,"upperquartile")
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
#maintain list of vertices in graph
vertexlist = sample_counts[:transcript_type]
##Measure of coexpression
#similarity_matrix=mutual_information(data)
#similarity_matrix = coexpression_measure(sample_data,"mutual_information")
similarity_matrix = coexpression_measure(sample_data,"PID")
## Adjacency matrix (using empricial distribution method atm)
threshold = 0.95
pre_adj_matrix = empirical_dist_adjacency(similarity_matrix,threshold)
#Trim nodes with degree zero
network_counts = sample_counts[vec(sum(pre_adj_matrix,dims=2).!=0),:]
network_data = sample_data[vec(sum(pre_adj_matrix,dims=2).!=0),:]

vertexlist = vertexlist[vec(sum(pre_adj_matrix,dims=2).!=0),:]
adj_matrix = copy(pre_adj_matrix)
adj_matrix = adj_matrix[:,vec(sum(pre_adj_matrix,dims=1).!=0)]
adj_matrix = adj_matrix[vec(sum(pre_adj_matrix,dims=2).!=0),:]
edgelist = edgelist_from_adj(adj_matrix)

#Network Analysis
degrees = sum(adj_matrix,dims=2)

##set up for distributed mode
#first clean to make sure there are no stray workers already around
using Distributed
rmprocs(workers())
#add workers equal to the number of available cpus	
addprocs(Threads.nthreads())
#addprocs(8)
@everywhere include(ENV["JULIA_PROJECT"]*"/src/GraphletCounting.jl")
@everywhere include(ENV["JULIA_PROJECT"]*"/src/NullModel.jl")
@time graphlet_counts = count_graphlets(vec(vertexlist),edgelist,4,"distributed")
#graphlet_concentrations = concentrate(graphlet_counts) 

@time motif_counts = find_motifs(edgelist,"hetero_rewire",100, typed = true, typelist = vec(vertexlist),plotfile="test.svg",graphlet_size = 4)
