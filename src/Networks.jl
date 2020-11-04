using LinearAlgebra
function adjacency(data,threshold)
	similarity_matrix=data
	similarity_matrix[diagind(similarity_matrix)].= 0
	similarity_matrix[broadcast(abs,similarity_matrix).<threshold].=0
	similarity_matrix[broadcast(abs,similarity_matrix).>threshold].=1
	similarity_matrix=BitArray(similarity_matrix)
	return similarity_matrix
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
