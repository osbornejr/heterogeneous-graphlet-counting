using DataStructures
#module GraphletCounting
#export Neighbours, mergecum, add3graphlets
#activate environment
##

##
#Functions
function Neighbours(edgelist::Array{Int,2})
	Neigh=Dict{Int,Vector{Int}}()
	for row in eachrow(edgelist)
       		vals=get!(Vector{Int},Neigh,row[1])
       		push!(vals,row[2])
       		vals=get!(Vector{Int},Neigh,row[2])
       		push!(vals,row[1])
       	end
	return Neigh 
end

function Neighbours(edgelist::Array{Pair{Int,Int},1})
	Neigh=Dict{Int,Vector{Int}}()
	for row in edgelist
       		vals=get!(Vector{Int},Neigh,row.first)
       		push!(vals,row.second)
		vals=get!(Vector{Int},Neigh,row.second)
       		push!(vals,row.first)
       	end
	return Neigh 
end

function Neighbours(edgelist::Array{Pair,1})
	Neigh=Dict{Int,Vector{Int}}()
	for row in edgelist
       		vals=get!(Vector{Int},Neigh,row.first)
       		push!(vals,row.second)
		vals=get!(Vector{Int},Neigh,row.second)
       		push!(vals,row.first)
       	end
	return Neigh 
end

function mergecum(d1::AbstractDict,d2::AbstractDict)
	return merge(+,d1,d2)
end

#improved version that should be faster? Using a default dict method
function add3graphlets(vertexlist::Array{String,1},nodelist::Array{Int,1},count_dict::DefaultDict{String,Int},i::Int,j::Int;graphlet_type::String)
	#type_col=2 ## may need to be more sophisticated at some point; assumes types are
		   ## given in the 2nd column of vertexlist. Alternatively, JUST give the type column.
	if length(nodelist)==0 
		return count_dict
	end
	delim="_"
	for (ind ,n) in enumerate(nodelist)
		count_dict[string(vertexlist[i],delim,vertexlist[j],delim,vertexlist[n],delim,graphlet_type)]+=1
	end 
	return count_dict
end

function graphlet_string(a::String,b::String,c::String,d::String,graphlet::String,delim::String)
	return string(a,delim,b,delim,c,delim,d,delim,graphlet)
end

function per_edge_counts(edge::Int,vertex_type_list::Array{String,1},edgelist::Union{Array{Pair{Int,Int},1},Array{Pair,1}},graphlet_size::Int,neighbourdict::Dict{Int,Vector{Int}},neighbourdictfunc::Function,ordered_vertices::Array{Int,1})
	count_dict = DefaultDict{String,Int}(0)
	h=edge	
#	# get nodes i and j for this edge
        i = edgelist[h].first
        j = edgelist[h].second
        #get neighbourhoods of i and j
        gamma_i = neighbourdict[i]	
        gamma_j = neighbourdict[j]
        
        #three node graphlets
        #Triangles: nodes connected to both i and j
        Tri = intersect(gamma_i,gamma_j)
        count_dict = add3graphlets(vertex_type_list,Tri,count_dict,i,j,graphlet_type="3-tri")
        # Paths with j at centre
	jPath = sort(setdiff(gamma_j,union(gamma_i,i)))
        count_dict = add3graphlets(vertex_type_list,jPath,count_dict,i,j,graphlet_type="3-path")
        # Paths with i at centre
	iPath = sort(setdiff(gamma_i,union(gamma_j,j)))
        count_dict = add3graphlets(vertex_type_list,iPath,count_dict,j,i,graphlet_type="3-path")
        if (graphlet_size==4)
        delim = "_"
	
	#four node graphlets
		for w in iPath
			for v in neighbourdict[w]
				if (v==i)

				elseif (!(v in gamma_i) & !(v in gamma_j))
						count_dict[graphlet_string(vertex_type_list[j],vertex_type_list[i],vertex_type_list[w],vertex_type_list[v],"4-path-edge-orbit",delim)]+=1
				elseif ((v in iPath) & (v < w))
						count_dict[graphlet_string(vertex_type_list[w],vertex_type_list[v],vertex_type_list[i],vertex_type_list[j],"4-tail-edge-orbit",delim)]+=1
 				end

			end
		end
						
		for w in jPath 
			for v in neighbourdict[w]
				if (v==j)
				#do nothing
				elseif (!(v in gamma_i) & !(v in gamma_j))
					count_dict[graphlet_string(vertex_type_list[i],vertex_type_list[j],vertex_type_list[w],vertex_type_list[v],"4-path-edge-orbit",delim)]+=1
				elseif ((v in jPath) & (v < w))
						count_dict[graphlet_string(vertex_type_list[w],vertex_type_list[v],vertex_type_list[j],vertex_type_list[i],"4-tail-edge-orbit",delim)]+=1
				elseif (v in iPath)
						count_dict[graphlet_string(vertex_type_list[v],vertex_type_list[i],vertex_type_list[j],vertex_type_list[w],"4-cycle",delim)]+=1
				end

			end
		end
	
		for w in Tri 
			for v  in neighbourdict[w]
			 	if (v==i|v==j)
				#do nothing
				elseif ((v in Tri) & (v < w)) 
					count_dict[graphlet_string(vertex_type_list[w],vertex_type_list[i],vertex_type_list[j],vertex_type_list[v],"4-clique",delim)]+=1
				## separating the processes here so that we can maintain the right type ordering 
				elseif (v in iPath) 
					count_dict[graphlet_string(vertex_type_list[v],vertex_type_list[i],vertex_type_list[w],vertex_type_list[j],"4-chord-edge-orbit",delim)]+=1
				elseif (v in jPath) 
					count_dict[graphlet_string(vertex_type_list[v],vertex_type_list[j],vertex_type_list[w],vertex_type_list[i],"4-chord-edge-orbit",delim)]+=1
				elseif (!(v in gamma_i) & !(v in gamma_j))
						count_dict[graphlet_string(vertex_type_list[i],vertex_type_list[j],vertex_type_list[w],vertex_type_list[v],"4-tail-tri-centre-orbit",delim)]+=1
				end				
			end
		end
		
		#Combinatorial methods
		# Remaining graphlets are found as per Rossi et al. (2019) using combinatorial relationships. These must be done per type 
		
		#first we store the types, as well as their occurences in each adjacent set (to save on recomputing)
		types = unique(vertex_type_list)
		iPathTypes = Array{Int}(undef,length(types))
		for(ind,t) in enumerate(types)
			iPathTypes[ind] = sum(vertex_type_list[iPath].==t)
		end
		jPathTypes = Array{Int}(undef,length(types))
		for(ind,t) in enumerate(types)
			jPathTypes[ind] = sum(vertex_type_list[jPath].==t)
		end
		TriTypes = Array{Int}(undef,length(types))
		for(ind,t) in enumerate(types)
			TriTypes[ind] = sum(vertex_type_list[Tri].==t)
		end

				##Now we loop per type combination
		for (inda,a) in enumerate(types)	
			for (indb,b) in enumerate(types[inda:end])
				if (a == b)
					##4-path centre orbit
					#order doesn't matter here
					count_dict[graphlet_string(a,vertex_type_list[i],vertex_type_list[j],b,"4-path-centre-orbit",delim)] += iPathTypes[inda]*jPathTypes[inda+indb-1] - count_dict[graphlet_string(a,vertex_type_list[i],vertex_type_list[j],b,"4-cycle",delim)]
					## 4-chord-centre-orbit
					count_dict[graphlet_string(a,vertex_type_list[i],vertex_type_list[j],b,"4-chord-centre-orbit",delim)] += Int(0.5*TriTypes[inda]*(TriTypes[inda]-1)) - count_dict[graphlet_string(a,vertex_type_list[i],vertex_type_list[j],b,"4-clique",delim)]

					##For the other two graphlets, we operate differently if i and j are of different types
					if (vertex_type_list[i]!=vertex_type_list[j])
						##4-star
						#To maintain type order here, we also have to separate. Note we enforce that the centre of the star is THIRD listed (in line with the 4-tail edge orbit layout)
						#i-centre
						count_dict[graphlet_string(a,b,vertex_type_list[i],vertex_type_list[j],"4-star",delim)] += Int(0.5*iPathTypes[inda]*(iPathTypes[inda]-1)) - 	count_dict[graphlet_string(a,b,vertex_type_list[i],vertex_type_list[j],"4-tail-edge-orbit",delim)]
						#j-centre
						count_dict[graphlet_string(a,b,vertex_type_list[j],vertex_type_list[i],"4-star",delim)] += Int(0.5*jPathTypes[inda]*(jPathTypes[inda]-1)) - 	count_dict[graphlet_string(a,b,vertex_type_list[j],vertex_type_list[i],"4-tail-edge-orbit",delim)]
						##4-tail tri-edge orbit
						#i-centre
						count_dict[graphlet_string(a,vertex_type_list[j],vertex_type_list[i],b,"4-tail-tri-edge-orbit",delim)] += TriTypes[inda]*iPathTypes[inda] - 	count_dict[graphlet_string(a,vertex_type_list[i],b,vertex_type_list[j],"4-chord-edge-orbit",delim)]
						#j-centre
						count_dict[graphlet_string(a,vertex_type_list[i],vertex_type_list[j],b,"4-tail-tri-edge-orbit",delim)] += TriTypes[inda]*jPathTypes[inda] - 	count_dict[graphlet_string(a,vertex_type_list[j],b,vertex_type_list[i],"4-chord-edge-orbit",delim)]

					else # when i and j are also of same type
						##4-star
						#To maintain type order here, we also have to separate. Note we enforce that the centre of the star is THIRD listed (in line with the 4-tail edge orbit layout)
						count_dict[graphlet_string(a,b,vertex_type_list[i],vertex_type_list[j],"4-star",delim)] += Int(0.5*iPathTypes[inda]*(iPathTypes[inda]-1))+Int(0.5*jPathTypes[inda]*(jPathTypes[inda]-1)) - 	count_dict[graphlet_string(a,b,vertex_type_list[i],vertex_type_list[j],"4-tail-edge-orbit",delim)]
						##4-tail tri-edge orbit
						count_dict[graphlet_string(a,vertex_type_list[i],vertex_type_list[j],b,"4-tail-tri-edge-orbit",delim)] += TriTypes[inda]*iPathTypes[inda] + TriTypes[inda]*jPathTypes[inda] - count_dict[graphlet_string(a,vertex_type_list[i],b,vertex_type_list[j],"4-chord-edge-orbit",delim)]
					end
				else # when a and b are of different types
					##4-path centre orbit
					#to maintain order here, we diverge from ROssi et al and calculate each orientation separately
					count_dict[graphlet_string(a,vertex_type_list[i],vertex_type_list[j],b,"4-path-centre-orbit",delim)] += iPathTypes[inda]*jPathTypes[inda+indb-1] - count_dict[graphlet_string(a,vertex_type_list[i],vertex_type_list[j],b,"4-cycle",delim)]
					count_dict[graphlet_string(b,vertex_type_list[i],vertex_type_list[j],a,"4-path-centre-orbit",delim)] += iPathTypes[inda+indb-1]*jPathTypes[inda] - count_dict[graphlet_string(b,vertex_type_list[i],vertex_type_list[j],a,"4-cycle",delim)]
					## 4-chord-centre-orbit
					count_dict[graphlet_string(a,vertex_type_list[i],vertex_type_list[j],b,"4-chord-centre-orbit",delim)] += TriTypes[inda]*TriTypes[inda+indb-1] - count_dict[graphlet_string(a,vertex_type_list[i],vertex_type_list[j],b,"4-clique",delim)] - count_dict[graphlet_string(b,vertex_type_list[i],vertex_type_list[j],a,"4-clique",delim)]

					##For the other two graphlets, we operate differently if i and j are of different types
					 if (vertex_type_list[i]!=vertex_type_list[j])
						##4-star
						#i-centre
						count_dict[graphlet_string(a,b,vertex_type_list[i],vertex_type_list[j],"4-star",delim)] += iPathTypes[inda]*iPathTypes[inda+indb-1] - count_dict[graphlet_string(a,b,vertex_type_list[i],vertex_type_list[j],"4-tail-edge-orbit",delim)] - count_dict[graphlet_string(b,a,vertex_type_list[i],vertex_type_list[j],"4-tail-edge-orbit",delim)]#unsure if both need to be subtracted here? TEST					
						#j-centre
						count_dict[graphlet_string(a,b,vertex_type_list[j],vertex_type_list[i],"4-star",delim)] += jPathTypes[inda]*jPathTypes[inda+indb-1] - count_dict[graphlet_string(a,b,vertex_type_list[j],vertex_type_list[i],"4-tail-edge-orbit",delim)] - count_dict[graphlet_string(b,a,vertex_type_list[j],vertex_type_list[i],"4-tail-edge-orbit",delim)]#unsure if both need to be subtracted here? TEST
						##4-tail tri-edge orbit
						#Note that here we split again! to make sure we get the right tail type (i.e. a or b)
						#i-centre, a tail
						count_dict[graphlet_string(b,vertex_type_list[j],vertex_type_list[i],a,"4-tail-tri-edge-orbit",delim)] += TriTypes[inda+indb-1]*iPathTypes[inda] - count_dict[graphlet_string(a,vertex_type_list[i],b,vertex_type_list[j],"4-chord-edge-orbit",delim)]
						#i-centre, b tail
						count_dict[graphlet_string(a,vertex_type_list[j],vertex_type_list[i],b,"4-tail-tri-edge-orbit",delim)] += TriTypes[inda]*iPathTypes[inda+indb-1] - count_dict[graphlet_string(b,vertex_type_list[i],a,vertex_type_list[j],"4-chord-edge-orbit",delim)]
						#j-centre, a tail
						count_dict[graphlet_string(b,vertex_type_list[i],vertex_type_list[j],a,"4-tail-tri-edge-orbit",delim)] += TriTypes[inda+indb-1]*jPathTypes[inda] - count_dict[graphlet_string(a,vertex_type_list[j],b,vertex_type_list[i],"4-chord-edge-orbit",delim)]
						#j-centre, b tail
						count_dict[graphlet_string(a,vertex_type_list[i],vertex_type_list[j],b,"4-tail-tri-edge-orbit",delim)] += TriTypes[inda]*jPathTypes[inda+indb-1] - count_dict[graphlet_string(b,vertex_type_list[j],a,vertex_type_list[i],"4-chord-edge-orbit",delim)]
					else
						##4-star
						count_dict[graphlet_string(a,b,vertex_type_list[i],vertex_type_list[j],"4-star",delim)] += iPathTypes[inda]*iPathTypes[inda+indb-1] +  jPathTypes[inda]*jPathTypes[inda+indb-1]  - count_dict[graphlet_string(a,b,vertex_type_list[i],vertex_type_list[j],"4-tail-edge-orbit",delim)] - count_dict[graphlet_string(b,a,vertex_type_list[i],vertex_type_list[j],"4-tail-edge-orbit",delim)]#unsure if both need to be subtracted here? TEST					
						##4-tail tri-edge orbit
						# We still split by and b tails
						#a tail
						count_dict[graphlet_string(b,vertex_type_list[i],vertex_type_list[j],a,"4-tail-tri-edge-orbit",delim)] += TriTypes[inda+indb-1]*iPathTypes[inda] + TriTypes[inda+indb-1]*jPathTypes[inda] - count_dict[graphlet_string(a,vertex_type_list[j],b,vertex_type_list[i],"4-chord-edge-orbit",delim)]
						#j-centre, b tail
						count_dict[graphlet_string(a,vertex_type_list[i],vertex_type_list[j],b,"4-tail-tri-edge-orbit",delim)] += TriTypes[inda]*iPathTypes[inda+indb-1] + TriTypes[inda]*jPathTypes[inda+indb-1] - count_dict[graphlet_string(b,vertex_type_list[j],a,vertex_type_list[i],"4-chord-edge-orbit",delim)]
					end
				end
	 		end
	 	end
 	end
	#combinatorial process currently adds 0 entries if no candidates exist. Not an issue per se, but makes readability on smaller graphs annoying. for now we tidy up at the end, but might be more efficient to do during combinatorial loop?
	for g in collect(keys(count_dict))[collect(values(count_dict)).==0]
		delete!(count_dict,g)
	end
	return count_dict 
end

function count_graphlets(vertex_type_list::Array{String,1},edgelist::Union{Array{Pair{Int,Int},1},Array{Pair,1}},graphlet_size::Int=3)


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
	
	#per edge process
	Threads.@threads for h in 1 :size(edgelist,1)
		Chi[h] = per_edge_counts(h,vertex_type_list,edgelist,graphlet_size,neighbourdict,neighbourdictfunc,ordered_vertices)        
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

function concentrate(graphlet_counts::Dict{String,Int})
	conc = Dict{String,Float64}()
	s = sum(collect(values(graphlet_counts)))
	for el in graphlet_counts
		conc[el.first] = el.second/s
	end
	return conc
end

function find_motifs(adjacency_matrix::AbstractArray,null_model::String,null_num::Int; typed::Bool=false, typelist::Array{String,1}=nothing,plotfile::String="DONOTPLOT")
	##Calculate null model counts
	edgelists = edgelists_null_model(adjacency_matrix,null_num,null_model,typelist)
	null_model_calc = null_model_counts(typelist,edgelists)
	null_model_df = null_model_dataframe(null_model_calc) 
	
	# calculate real network counts
	graphlet_counts = count_graphlets(vertexlist[:,2],edgelist,4)
	
	#Statistical significance
	zscores = Dict{String,Float64}()
	for g in collect(keys(graphlet_counts))
		## get score that each null model graph got for corresponding graphlet. If it didn't appear at all in a random graph, we must add the 0 here as well
		rand_vals = filter(:graphlet=>x->x==g,null_model_df)[!,:value]
		rand_vals = vcat(rand_vals,zeros(Int,null_num-length(rand_vals)))
		m = mean(rand_vals)
		sd = std(rand_vals)
		Z = (graphlet_counts[g]-m)/sd
		zscores[g] = Z
	end
	
	if (plotfile !="DONOTPLOT")	
		#convert to concentrations for plotting
		null_model_conc = null_model_concentrations(null_model_calc)
		null_model_df = null_model_dataframe(null_model_conc) 
		graphlet_concentrations = concentrate(graphlet_counts) 
		graphlet_df = DataFrame(graphlet = broadcast(first,collect(graphlet_concentrations)),value = broadcast(last,collect(graphlet_concentrations)))
		p = plot(layer(null_model_df,x=:graphlet,y=:value,Geom.boxplot(suppress_outliers = true),color=:graphlet),layer(graphlet_df,x = :graphlet,y = :value, Geom.point,color=["count in graph"]),Guide.xticks(label=false));
		draw(SVG(plotfile,12inch,6inch),p)
		print("Stat plot saved to $plotfile.")
	end
	#null_model_sum = reduce(mergecum,null_model_calc)
	return zscores
end	
