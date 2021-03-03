using DataFrames, DataStructures,ProgressMeter,Distributed

#aggregator function
function t2(d1,d2)
	append!(d1,d2)
	d1
		
end

function edgelists_null_model(edgelist::Union{Array{Pair{Int,Int},1},Array{Pair,1}},n::Int,method::String,typelist::Array{String,1},graphlet_size::Int=3)
	#degrees = sum(adj_matrix,dims=2)
	#edgelist = edgelist_from_adj(adj_matrix) 
	edgelists = Array{Union{Array{Pair{Int,Int},1},Array{Pair,1}}}(undef,n)
	##using default switching value of 100*m as per Milo et al TODO user can set this
	switching_factor =100
	switching_val = switching_factor*length(edgelist)
	if (method == "erdos_renyi")
		for i in 1:n
       			edgelists[i] = edgelist_from_adj(LightGraphs.LinAlg.adj_matrix(erdos_renyi(size(degrees,1),Int(sum(degrees)/2))))
       		end
	end
	if(method == "switching")
		for i in 1:n
			print("Switching network $i...\n")
			edgelists[i] = edge_switch(edgelist,switching_val)  
	 	end
	end 
	if(method == "rewire")
		for i in 1:n	
			print("Rewiring network $i...\n")
			edgelists[i] = edgelist_from_adj(rewire(adj_matrix,switching_val))
			print("Successfully rewired edges.\n")	
		end
	end
	if(method == "hetero_rewire")
	       	@info "Rewiring $n graphs..."
		#Distributed method
		#edgelists = @showprogress @distributed (t2) for i in 1:n	
		#	[hetero_rewire(edgelist,switching_factor,typelist,graphlet_size)]
		#end

		for i in 1:n
			print("Switching network $i...\n")
			hetero_rewire(edgelist,switching_factor,typelist,graphlet_size)
	 	end
	end

	return edgelists
end


function null_model_counts(typelist::Array{String,1},edgelists::Array{Array{Pair,1},1}) 
	null_num = size(edgelists,1)
	null_model = Array{Dict{String,Int},1}(undef,null_num)
	@info "Counting null graphs..."
	for (i,el) in enumerate(edgelists)
		@info "For null graph $i:" 
		null_model[i] = count_graphlets(typelist,el,4,"distributed")
	end
	return null_model
end

function null_model_concentrations(null_model_counts::Array{Dict{String,Int},1})
##get concentrations from a set of calculated null_model_counts.
	null_num = size(null_model_counts,1)
	null_model = Array{Dict{String,Float64},1}(undef,null_num)
	for (i,mod) in enumerate(null_model_counts)
		null_dict = Dict{String,Float64}()
		s = sum(collect(values(mod)))
		for ent in mod
			null_dict[ent.first] = ent.second/s
		end
		null_model[i] = null_dict
	end

	return null_model
end

function null_model_dataframe(null_model_dicts::Array{Dict{String,T},1} where T<:Real)
#use this function to prepare null model for plot function, as a long format dataframe
	null_collection = union(null_model_dicts...)
	df = DataFrame(graphlet = broadcast(first,null_collection),value = broadcast(last,null_collection))
	return df
end

##allows a function to be broadcast and splat each argument
splat(f) = args->f(args...)

function switch_edges(edgelist::Union{Array{Pair{Int,Int},1},Array{Pair,1}},switches::Int)
##This works (seemingly) but it does not scale at all well to bigger switch requirements, which will be necessary on larger networks if we want num_of_switches = 100*m 
	#convert edgelist to a two column array of ints (i.e. mutable)
	el = hcat(first.(edgelist),last.(edgelist))
	#function that checks if new edges already exist (to prevent multiedges)
	check(edges,ihead,jtail) = sum(((view(edges,1:size(edges,1),1).==ihead)+(view(edges,1:size(edges,1),2).==jtail)).>1)==0 && sum(((view(edges,1:size(edges,1),1).==jtail)+(view(edges,1:size(edges,1),2).==ihead)).>1)==0;
	suc = 0
	for s in 1:switches
		#choose 2 edges at random
		i = rand(1:length(edgelist))
		j = rand(1:length(edgelist))
		if(i!=j)
			##if edge ends aren't shared-- preventing self-loops
			if (el[i,1]!=el[j,2] && el[j,1]!=el[i,2])
				
				#if both new edges do not already exist in el (prevent multiloops)
				if(check(el,el[i,1],el[j,2]) &&check(el,el[j,1],el[i,2])) 
				   	#switch ends
					temp = copy(el[i,2])
					el[i,2] = el[j,2]
					el[j,2] = temp
					suc += 1
				end
			end	
		end
	end
	#convert switched edges back to immutable Pair type
	print("Successfully switched $suc edges in network.\n")
	return output = splat(Pair).(eachrow(el))
end


function configuration_model(edgelistt::Union{Array{Pair{Int,Int},1},Array{Pair,1}})	
	#this is a work in progress, and probably not usable! Impossible to find a network without a self loops, let alone also without multiedges! Note that there are \Prod_i k_i! possible networks in the space of all graphs with the same degree distribution as G-- for the 10% network, this corresponds to a search space of ~10^430000! 
	count = 0
       	loops = 1
	stubs = hcat(first.(edgelist),last.(edgelist))
       	while (loops > 0)
       		shuffle!(stubs)
       		count += 1
       		loops = sum(stubs[:,1].==stubs[:,2]) 
       		if (loops < 30)
       			@info "Graph $count had $loops self-loops"
       		end
       		if (mod(count,10000)==0)
       			@info "Checked $count graphs"
       		end
       	end
end

function edge_switch(edgelist::Union{Array{Pair{Int,Int},1},Array{Pair,1}},switches::Int,hetero::Bool = false)
	#convert edgelist to dictionary form
	el = Dict{Pair,Bool}()
	for e in edgelist
		el[e] = 1
	end
	diffo = 0
	for s in 1:switches
		i = first(rand(el))
		j = first(rand(el))
		if(i!=j)
			if(hetero == true)
				##in this case self-loops aren't a concern; by definition the different typed ends will be unique
				
				##MULTI-EDGES: check putative new edges aren't already in dict.
				if(!haskey(el,Pair(first(i),last(j))) && !haskey(el,Pair(first(j),last(i))) && !haskey(el,Pair(last(j),first(i))) && !haskey(el,Pair(last(i),first(j))))
					delete!(el,i)
					delete!(el,j)
					## just switch ahead in the only way possible (a la a pseudo directed network)
					el[Pair(first(i),last(j))] = 1
					el[Pair(first(j),last(i))] = 1
				end
			else  #hetero = false i.e. homogeneous case
	
				#to maintain unbiased sampling, we need to switch either fronts to ends or fronts to fronts (and ends to ends) with equal probability
				if (rand()>0.5)	#starts to ends
					##SELF-LOOPS: check ends aren't shared in i and j
			 		if(first(i)!=last(j) && first(j)!=last(i))
					##MULTI-EDGES: check putative new edges aren't already in dict
					 	if(!haskey(el,Pair(first(i),last(j))) && !haskey(el,Pair(first(j),last(i))) && !haskey(el,Pair(last(i),first(j))) && !haskey(el,Pair(last(j),first(i))))
						 	delete!(el,i)
					 		delete!(el,j)
							el[Pair(first(i),last(j))] = 1
				 			el[Pair(first(j),last(i))] = 1
						end
					end
				else #starts to starts, ends to ends
					##SELF-LOOPS: check starts and ends aren't shared in i and j
			 		if(first(i)!=first(j) && last(j)!=last(i))
					##MULTI-EDGES: check putative new edges aren't already in dict
					 	if(!haskey(el,Pair(first(i),first(j))) && !haskey(el,Pair(last(j),last(i))) && !haskey(el,Pair(first(j),first(i))) && !haskey(el,Pair(last(i),last(j))))
						 	delete!(el,i)
					 		delete!(el,j)
							el[Pair(first(i),first(j))] = 1
				 			el[Pair(last(i),last(j))] = 1
						end
				 	end
				end

			end
		end
	end
	return collect(keys(el))
end

### function to quickly implement r-igraph method in julia (as it is ~10 times faster than best method I could come up with in pure julia at present). 
#function rewire(adj::AbstractArray,switching_val::Int)
#	
#	@rput adj
#	@rput switching_val
#	print("Rewiring network...\n")
#	R"""
#	library(igraph)
#	g = graph.adjacency(adj,mode="undirected",diag=F)
#	ng = rewire(g,keeping_degseq(loops=FALSE,niter=switching_val))
#	a = as.matrix(get.adjacency(ng));
#	"""
#	adjacency = BitArray(@rget a)
#end	

#uses the above rewire function to only switch edges that share the same typed end points. In this case, a switching factor is provided such that switching_factor*m = switching_val TODO maybe add this approach to all rewiring/switching algorithms? default 100
function hetero_rewire(edgelist::Union{Array{Pair{Int,Int},1},Array{Pair,1}},switching_factor::Int,typelist::Array{String,1},graphlet_size::Int=3)
	##first subset the edgelist based on the unique endpoint pairs for edges 
	edgetypes = splat(Pair).(eachrow(hcat(map(x-> typelist[x],first.(edgelist)),map(x-> typelist[x],last.(edgelist)))))
	t = length(unique(edgetypes))	
	@info "Found the unique edgetypes: $(unique(edgetypes))"
	#now we form a subnetwork for each "edgetype", and rewire that. The results are 
	#cumulatively stored in an intitially empty edgelist 
	adjusted_edges = Array{Pair,1}()
	for type in unique(edgetypes)
		s_edges = edgelist[edgetypes.==type]
		## to maintain uniformity, the switching algorithm should be implemented differently depending on whether the edgetype is homogeneous or heterogeneous
		##homogeneous case
		if (first(type)==last(type))
                         s_adjust = edge_switch(s_edges,switching_factor*length(s_edges))
		 	 adjusted_edges = vcat(adjusted_edges,s_adjust)
		else
		##heterogeneous case
                         s_adjust = edge_switch(s_edges,switching_factor*length(s_edges),true)
		 	 adjusted_edges = vcat(adjusted_edges,s_adjust)
		end
	end
	if (graphlet_size == 4)
		@info "Matching random graph to target 3-node graphlet vector..." 
		#this involves calculating the 3 node subgraphs at every step (laborious! and probably not practical at present).
		randomised = DefaultDict(0,count_graphlets(typelist,adjusted_edges)[1])
		target = DefaultDict(0,count_graphlets(typelist,edgelist)[1])
		##define energy on difference between adjusted and target counts
		E = sum([abs(randomised[g]-target[g])/(randomised[g]+target[g]) for g in collect(keys(target))])
		Temp = 10
		edgetypes = splat(Pair).(eachrow(hcat(map(x-> typelist[x],first.(adjusted_edges)),map(x-> typelist[x],last.(adjusted_edges)))))
		@info "Need zero energy:"
		while (E>0.05)
			trial = Array{Pair,1}()
			## try an edge switch (we choose an edgetype and only switch in that subgraph)
			chosen_type = rand(unique(edgetypes))
			for type in unique(edgetypes)
				s_edges = adjusted_edges[edgetypes.==type]
				## to maintain uniformity, the switching algorithm should be implemented differently depending on whether the edgetype is homogeneous or heterogeneous
				##homogeneous case
				if (type == chosen_type)
					if (first(type)==last(type))
        		                	 s_adjust = edge_switch(s_edges,1)
				 	 	trial = vcat(trial,s_adjust)
					else
				##heterogeneous case
        		                	s_adjust = edge_switch(s_edges,1,true)
				 	 	trial = vcat(trial,s_adjust)
					end
				else
					trial = vcat(trial,s_edges)
				end
			end
			#calculate new 3-node counts
			randomised = DefaultDict(0,count_graphlets(typelist,trial)[1])
			#E trial
			E_trial = sum([abs(randomised[g]-target[g])/(randomised[g]+target[g]) for g in collect(keys(target))])
			## update if change is accepted
			crit = rand()
			if (E_trial<E)
				adjusted_edges = trial
				edgetypes = splat(Pair).(eachrow(hcat(map(x-> typelist[x],first.(adjusted_edges)),map(x-> typelist[x],last.(adjusted_edges)))))
				E = E_trial
				@info "Energy is $E, temperature is $Temp"
				#cool slightly
				Temp = Temp*0.99
			elseif (crit<exp(-(1+E_trial-E)/Temp))
				prob = exp(-(1+E_trial-E)/Temp)
				adjusted_edges = trial
				edgetypes = splat(Pair).(eachrow(hcat(map(x-> typelist[x],first.(adjusted_edges)),map(x-> typelist[x],last.(adjusted_edges)))))
				E = E_trial
				@info "Increased energy is $E, temperature is $Temp. Rand was $crit, probability was $prob."
				#cool slightly
				Temp = Temp*0.99


			end
			##reheat if solidified too much
			if (Temp<1e-300)
				Temp = 1e-100
			end
		end
	end
	return adjusted_edges
end

#THIS VERSION IS NOT WORKING AT PRESENT, but should be a faster option if it can be got to work using the igraph based rewire function.
function hetero_rewire(adj_matrix::AbstractArray,switching_factor::Int,typelist::Array{String,1})
	##first subset the edgelist based on the unique endpoint pairs for edges 
	print("Partitioning network based on edge endpoint types...\n")
	edgelist = edgelist_from_adj(adj_matrix)
	edgetypes = splat(Pair).(eachrow(hcat(map(x-> typelist[x],first.(edgelist)),map(x-> typelist[x],last.(edgelist)))))
	t = length(unique(edgetypes))	
	print("Found $t unique edgetpyes to partition on:\n")
	print(unique(edgetypes))
	#now we form a subnetwork for each "edgetype", and rewire that. The results are 
	#cumulatively stored in an intitially empty copy of the input adjacency matrix
	adjusted = copy(adj_matrix)
	for type in unique(edgetypes)
		s_edges = edgelist[edgetypes.==type]
		s_nodes = sort(unique(vcat(first.(s_edges),last.(s_edges))))
		s_adj_matrix = adj_matrix[s_nodes,s_nodes]
		s_adjusted = rewire(s_adj_matrix,switching_factor*Int(sum(s_adj_matrix)/2))
		adjusted[s_nodes,s_nodes] = s_adjusted
	end
	return adjusted
end

