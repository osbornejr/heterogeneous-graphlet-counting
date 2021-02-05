using MPI

function count_graphlets_mpi(vertex_type_list::Array{String,1},edgelist::Union{Array{Pair{Int,Int},1},Array{Pair,1}},graphlet_size::Int=3)


	##INPUTS TO PER EDGE FUNCTION
	#get neighbourhood for each vertex in advance (rather than calling per-edge)
	neighbourdict=Neighbours(edgelist)
	# set up function to apply dict to arrays (maybe a better way exists, but this works?)
	neighbourdictfunc(x::Int) = neighbourdict[x]
	#vertex set derived from edgelist (allows for any arbritary vertex labelling)
	##(important: not an ordered list that matches with vertexlist types!!) 
	ordered_vertices = unique(vcat(first.(edgelist),last.(edgelist)))		
	
	
	
	#preallocate array to store each edge's graphlet dictionary 
	Chi=Array{Dict}(undef,size(edgelist,1));
	
	#per edge process: attempting to use MPI here.
	
	MPI.Init()
	comm = MPI.COMM_WORLD
	rank = MPI.Comm.rank(comm)
	nprocs = MPI.Comm.size(comm)
	#crude method to partition sets of edges to each available process. 
	count = Int(floor(size(edgelist)/nprocs))
	if rank==(nprocs-1) ##last process can do excess edges
		start = rank*count
		stop = size(edgelist) 
	else
		start = rank*count
		stop = start + count
	end
	
	for h in 1:size(edgelist,1) 
		count_dict = per_edge_counts(h,vertex_type_list,edgelist,graphlet_size,neighbourdict,neighbourdictfunc,ordered_vertices)        
		Chi[h] = count_dict		
	end
	#total counts for each graphlet
	total_counts = reduce(mergecum,Chi)
	
	#reorder names to merge orbits
	graphlet_names = (split.(collect(keys(total_counts)),"_"))
	for el in 1:size(graphlet_names,1)
		#for 3 graphlets:
		if(length(graphlet_names[el])==4)
		   
			#do not reorder for x-y-x paths (different orbit to other 3-paths) 
			if !(graphlet_names[el][1]!=graphlet_names[el][2] && graphlet_names[el][1]==graphlet_names[el][3] && graphlet_names[el][4]=="3-path")
	          		graphlet_names[el][1:3]=sort(graphlet_names[el][1:3])
 	        	end
		end
		## for 4 graphlets
		if(length(graphlet_names[el])==5)
			#clean off orbit listing
			graphlet_names[el][5] = replace(graphlet_names[el][5],Pair("-edge-orbit",""))
			graphlet_names[el][5] = replace(graphlet_names[el][5],Pair("-centre-orbit",""))
			graphlet_names[el][5] = replace(graphlet_names[el][5],Pair("-tri",""))
			
			#paths (maintain and order centre edge, moving others accordingly)
			if (graphlet_names[el][5] == "4-path")
				#we iwant to switch if interior needs switching: 
				if (graphlet_names[el][[2,3]][1] != sort(graphlet_names[el][[2,3]])[1])
					graphlet_names[el][[2,3]] = sort(graphlet_names[el][[2,3]])
					graphlet_names[el][[1,4]] = graphlet_names[el][[4,1]]
				end
				#or if interior is the same, we sort outer types:
				if (graphlet_names[el][[2]] == graphlet_names[el][[3]])
					graphlet_names[el][[1,4]] = sort(graphlet_names[el][[1,4]]) 
				end
			#stars (maintain star centre (3rd entry), order others)
			elseif (graphlet_names[el][5] == "4-star")
				graphlet_names[el][[1,2,4]] = sort(graphlet_names[el][[1,2,4]])
			#tails (maintain and order edge not connected to tail)
			elseif(graphlet_names[el][5] == "4-tail")
				graphlet_names[el][[1,2]] = sort(graphlet_names[el][[1,2]])
			#chords (maintain and order centre edge in middle of name)
			elseif (graphlet_names[el][5] == "4-chord")
				graphlet_names[el][[2,3]] = sort(graphlet_names[el][[2,3]])
				graphlet_names[el][[1,4]] = sort(graphlet_names[el][[1,4]])
			#cycles and cliques (just order everything)
			else 
				graphlet_names[el][1:4] = sort(graphlet_names[el][1:4])
			end
		end

 	end
	#merge counts for each graphlet type
	graphlet_counts = Dict{String,Int}()
	orbits = join.(graphlet_names,"_")
	count_values = collect(values(total_counts))
	for orb in 1:size(orbits,1)
		graphlet_counts[orbits[orb]] = get(graphlet_counts,orbits[orb],0)+count_values[orb]	
	end
	#divide each graphlet count by number of edges in graphlet
	for g in collect(keys(graphlet_counts))
		if (occursin("3-tri",g))
			graphlet_counts[g] = div(graphlet_counts[g],3)
		end
		if (occursin("3-path",g))
			graphlet_counts[g] = div(graphlet_counts[g],2)
		end
		if (occursin("4-path",g))
			graphlet_counts[g] = div(graphlet_counts[g],3)
		end
		if (occursin("4-star",g))
			graphlet_counts[g] = div(graphlet_counts[g],3)
		end
		if (occursin("4-cycle",g))
			graphlet_counts[g] = div(graphlet_counts[g],4)
		end
		if (occursin("4-tail",g))
			graphlet_counts[g] = div(graphlet_counts[g],4)
		end
		if (occursin("4-chord",g))
			graphlet_counts[g] = div(graphlet_counts[g],5)
		end
		if (occursin("4-clique",g))
			graphlet_counts[g] = div(graphlet_counts[g],6)
		end
	end
	return graphlet_counts
end

MPI.Finalize()
