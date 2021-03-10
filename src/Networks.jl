using LinearAlgebra
function adjacency(data,threshold)
	sim_matrix = copy(data)
	sim_matrix[diagind(sim_matrix)].= 0
	sim_matrix[broadcast(abs,sim_matrix).<threshold].=0
	sim_matrix[broadcast(abs,sim_matrix).>threshold].=1
	sim_matrix=BitArray(sim_matrix)
	return sim_matrix
end

function edgelist_from_adj(adjacency_matrix)
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
