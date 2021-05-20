using LinearAlgebra
function adjacency(data::AbstractArray,threshold::Float64)
	sim_matrix = copy(data)
	sim_matrix[diagind(sim_matrix)].= 0
	sim_matrix[broadcast(abs,sim_matrix).<threshold].=0
	sim_matrix[broadcast(abs,sim_matrix).>threshold].=1
	sim_matrix=BitArray(sim_matrix)
	return sim_matrix
end
function empirical_dist_adjacency(sim_matrix::AbstractArray,prob::Float64)
	##run each row of similarity matrix through its own empirical cdf to determine probability 
	#dists = hcat(map.(ecdf.(eachrow(sim_matrix)),eachrow(sim_matrix))...)
	##alternative method to ignore 0 values in ecdfs...
	dists = hcat(map.(ecdf.(filter.(x->x!=0,eachrow(sim_matrix))),eachrow(sim_matrix))...)
	##find which are significant for each row and column, then combine
	adj = (dists.>prob).*(dists'.>prob)
	adj[diagind(adj)].= 0
	return adj
end

function top_adjacency(sim_matrix::AbstractArray,top_scores::Int)
	sim_matrix[diagind(sim_matrix)].= 0
	##sort each row of similarity matrix
	sorted = map(sort,eachrow(broadcast(abs,sim_matrix)))
	## get lowest included score out of top scores for each row
	lims = getindex.(sorted,Ref(size(sim_matrix,1)-top_scores))
	cuts = broadcast(abs,sim_matrix).>lims
	##check against transpose to see which entries agree for BOTH transcripts
	adj = cuts.*cuts'
	return adj
end

function edgelist_from_adj(adjacency_matrix::AbstractArray)
	edgelist=Array{Pair}(undef,sum(UpperTriangular(adjacency_matrix)))
	count=0
	for (i,row) in enumerate(eachrow(UpperTriangular(adjacency_matrix)))
       		for j in 1:size(row,1)
			if (row[j]==1)
				count=count+1
				edgelist[count]=Pair(i,j)
       			end
       		end
       end
	return edgelist
end
