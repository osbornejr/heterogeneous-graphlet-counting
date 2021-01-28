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



function add3graphlets(vertexlist::Array{String,1},nodelist::Array{Int,1},count_dict::Dict{String,Int},i::Int,j::Int;graphlet_type::String)
	#type_col=2 ## may need to be more sophisticated at some point; assumes types are
		   ## given in the 2nd column of vertexlist. Alternatively, JUST give the type column.
	if length(nodelist)==0 
		return count_dict
	end
	graphlet_list=fill("",length(nodelist))
	delim="_"
	for (ind ,n) in enumerate(nodelist)
		graphlet_list[ind]=string(vertexlist[i],delim,vertexlist[j],delim,vertexlist[n],delim,graphlet_type)
	end 
	#list=vertexlist[i].*"_".*vertexlist[j].*"_".*vertexlist[nodelist].*"_".*graphlet_type	
	x=Dict{String,Integer}([(i,count(x->x==i,graphlet_list)) for i in unique(graphlet_list)]) 
	merge!(+,count_dict,x)
	return count_dict
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

function add4graphlets(typelist::Array{String,1},nodelist1::Array{Int,1},nodelist2::Array{Array{Int,1},1},count_dict::Dict{String,Int},i::Int,j::Int;graphlet_type::String,order::Array{Int,1}=[1,2,3,4])

	#check to see if candidate lists are empty
	if length(nodelist1)==0 
		return count_dict
	elseif sum(length.(nodelist2))==0
		return count_dict
	else ##we do have some graphlets to add!
		delim = "_"
		for (ind,n) in enumerate(nodelist1)
			#check if specific list is empty
			if length(nodelist2[ind])==0
				x = Dict{String,Int}()
			else
				list = fill("",sum(length(nodelist2[ind])))
				if (graphlet_type=="4-path-edge-orbit")
					for (indd,m) in enumerate(nodelist2[ind])
						#ordering i_j_jneigh_j_neighneigh
						list[indd] = string(typelist[i],delim,typelist[j],delim,typelist[n],delim,typelist[m],delim,graphlet_type)
					end
				end
				if (graphlet_type == "4-path-centre-orbit")
					for (indd,m) in enumerate(nodelist2[ind])
						#ordering neigh_i_i_j_jneigh
						list[indd] = string(typelist[m],delim,typelist[i],delim,typelist[j],delim,typelist[n],delim,graphlet_type)
					end
				end
                                if (graphlet_type == "4-tail-edge-orbit" )
					for (indd,m) in enumerate(nodelist2[ind])
						#ordering neighi_neighi_i_j
                                                list[indd] = string(typelist[n],delim,typelist[m],delim,typelist[i],delim,typelist[j],delim,graphlet_type)
					end
				end
                                if (graphlet_type == "4-star")
					for (indd,m) in enumerate(nodelist2[ind])
						#ordering j_neighi_neighi_i
                                                list[indd] = string(typelist[j],delim,typelist[n],delim,typelist[m],delim,typelist[i],delim,graphlet_type)
					end
                                        
                                end
                                if (graphlet_type == "4-tail-tri-edge-orbit" )
					for (indd,m) in enumerate(nodelist2[ind])
						#ordering k_j_i_neigh_i
						list[indd] = string(typelist[m],delim,typelist[j],delim,typelist[i],delim,typelist[n],delim,graphlet_type)
					end
                                        
                                end
                                if (graphlet_type == "4-tail-tri-centre-orbit" )
					for (indd,m) in enumerate(nodelist2[ind])
						#ordering i_j_k_neigh_k
						list[indd] = string(typelist[i],delim,typelist[j],delim,typelist[n],delim,typelist[m],delim,graphlet_type)
					end
                                        
                                end
                                if (graphlet_type == "4-cycle")
					for (indd,m) in enumerate(nodelist2[ind])
						#ordering neighi_i_j_neighj
                                                list[indd] = string(typelist[n],delim,typelist[i],delim,typelist[j],delim,typelist[m],delim,graphlet_type)
					end
                                        
                                end
                                if (graphlet_type == "4-chord-edge-orbit")
					for (indd,m) in enumerate(nodelist2[ind])
						#ordering neighik_i_k_j
                                                list[indd] = string(typelist[m],delim,typelist[i],delim,typelist[n],delim,typelist[j],delim,graphlet_type)
					end
                                        
                                end
                                if (graphlet_type == "4-chord-centre-orbit")
					for (indd,m) in enumerate(nodelist2[ind])
						#ordering neighij_i_j_k
                                                list[indd] = string(typelist[m],delim,typelist[i],delim,typelist[j],delim,typelist[n],delim,graphlet_type)
					end
                                        
                                end
                                if (graphlet_type == "4-clique")
					for (indd,m) in enumerate(nodelist2[ind])
						#ordering k_i_j_neigh_kij
                                                list[indd] = string(typelist[n],delim,typelist[i],delim,typelist[j],delim,typelist[m],delim,graphlet_type)
					end
                                        
                                end
				x = Dict{String,Int}([(p,count(x->x==p,list)) for p in unique(list)])
			end
		merge!(+,count_dict,x)
		end
	end
	return count_dict
end


#improved version that should be faster? Using a default dict method
function add4graphlets(typelist::Array{String,1},nodelist1::Array{Int,1},nodelist2::Array{Array{Int,1},1},count_dict::DefaultDict{String,Int},i::Int,j::Int;graphlet_type::String,order::Array{Int,1}=[1,2,3,4])

	#check to see if candidate lists are empty
	if length(nodelist1)==0 
		return count_dict
	elseif sum(length.(nodelist2))==0
		return count_dict
	else ##we do have some graphlets to add!
		delim = "_"
		for (ind,n) in enumerate(nodelist1)
			#check if specific list is empty
			if length(nodelist2[ind])==0
			#do nothing 
			else
				if (graphlet_type=="4-path-edge-orbit")
					for (indd,m) in enumerate(nodelist2[ind])
						#ordering i_j_jneigh_j_neighneigh
						count_dict[string(typelist[i],delim,typelist[j],delim,typelist[n],delim,typelist[m],delim,graphlet_type)]+=1
					end 
				end
				if (graphlet_type == "4-path-centre-orbit")
					for (indd,m) in enumerate(nodelist2[ind])
						#ordering neigh_i_i_j_jneigh
						count_dict[string(typelist[m],delim,typelist[i],delim,typelist[j],delim,typelist[n],delim,graphlet_type)]+=1
					end
				end
                                if (graphlet_type == "4-tail-edge-orbit" )
					for (indd,m) in enumerate(nodelist2[ind])
						#ordering neighi_neighi_i_j
						count_dict[string(typelist[n],delim,typelist[m],delim,typelist[i],delim,typelist[j],delim,graphlet_type)]+=1
					end
				end
                                if (graphlet_type == "4-star")
					for (indd,m) in enumerate(nodelist2[ind])
						#ordering j_neighi_neighi_i
						count_dict[string(typelist[j],delim,typelist[n],delim,typelist[m],delim,typelist[i],delim,graphlet_type)]+=1
					end
                                        
                                end
                                if (graphlet_type == "4-tail-tri-edge-orbit" )
					for (indd,m) in enumerate(nodelist2[ind])
						#ordering k_j_i_neigh_i
						count_dict[string(typelist[m],delim,typelist[j],delim,typelist[i],delim,typelist[n],delim,graphlet_type)]+=1
					end
                                        
                                end
                                if (graphlet_type == "4-tail-tri-centre-orbit" )
					for (indd,m) in enumerate(nodelist2[ind])
						#ordering i_j_k_neigh_k
						count_dict[string(typelist[i],delim,typelist[j],delim,typelist[n],delim,typelist[m],delim,graphlet_type)]+=1
					end
                                        
                                end
                                if (graphlet_type == "4-cycle")
					for (indd,m) in enumerate(nodelist2[ind])
						#ordering neighi_i_j_neighj
						count_dict[string(typelist[n],delim,typelist[i],delim,typelist[j],delim,typelist[m],delim,graphlet_type)]+=1
					end
                                        
                                end
                                if (graphlet_type == "4-chord-edge-orbit")
					for (indd,m) in enumerate(nodelist2[ind])
						#ordering neighik_i_k_j
						count_dict[string(typelist[m],delim,typelist[i],delim,typelist[n],delim,typelist[j],delim,graphlet_type)]+=1
					end
                                        
                                end
                                if (graphlet_type == "4-chord-centre-orbit")
					for (indd,m) in enumerate(nodelist2[ind])
						#ordering neighij_i_j_k
						count_dict[string(typelist[m],delim,typelist[i],delim,typelist[j],delim,typelist[n],delim,graphlet_type)]+=1
					end
                                        
                                end
                                if (graphlet_type == "4-clique")
					for (indd,m) in enumerate(nodelist2[ind])
						#ordering k_i_j_neigh_kij
						count_dict[string(typelist[n],delim,typelist[i],delim,typelist[j],delim,typelist[m],delim,graphlet_type)]+=1
					end
                                        
                                end
			end
		end
	end
	return count_dict
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
        jPath = setdiff(gamma_j,union(gamma_i,i))
        count_dict = add3graphlets(vertex_type_list,jPath,count_dict,i,j,graphlet_type="3-path")
        # Paths with i at centre
        iPath = setdiff(gamma_i,union(gamma_j,j))
        count_dict = add3graphlets(vertex_type_list,iPath,count_dict,j,i,graphlet_type="3-path")
        if (graphlet_size==4)
        delim = "_"
		
		for w in iPath
			for v in neighbourdict[w]
				if (v==i|v==j)

				elseif (!(v in gamma_i) & !(v in gamma_j))
						count_dict[string(vertex_type_list[j],delim,vertex_type_list[i],delim,vertex_type_list[w],delim,vertex_type_list[v],delim,"4-path-edge-orbit")]+=1
				elseif ((v in iPath) & (v < w))
						count_dict[string(vertex_type_list[w],delim,vertex_type_list[v],delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-tail-edge-orbit")]+=1
				end

			end
		end
						
		for w in jPath
			for v in neighbourdict[w]
				if (v==i|v==j)
				#do nothing
				elseif (!(v in gamma_i) & !(v in gamma_j))
						count_dict[string(vertex_type_list[i],delim,vertex_type_list[j],delim,vertex_type_list[w],delim,vertex_type_list[v],delim,"4-path-edge-orbit")]+=1
				elseif ((v in jPath) & (v < w))
						count_dict[string(vertex_type_list[w],delim,vertex_type_list[v],delim,vertex_type_list[j],delim,vertex_type_list[i],delim,"4-tail-edge-orbit")]+=1
				elseif (v in iPath)
						count_dict[string(vertex_type_list[w],delim,vertex_type_list[i],delim,vertex_type_list[j],delim,vertex_type_list[v],delim,"4-cycle")]+=1
				end

			end
		end
	
		for w in Tri 
			for v  in neighbourdict[w]
			 	if (v==i|v==j)
				#do nothing
				elseif ((v in Tri) & (v < w)) 
					count_dict[string(vertex_type_list[w],delim,vertex_type_list[i],delim,vertex_type_list[j],delim,vertex_type_list[v],delim,"4-clique")]+=1
				## separating the processes here so that we can maintain the right type ordering 
				elseif (v in iPath) 
					count_dict[string(vertex_type_list[v],delim,vertex_type_list[i],delim,vertex_type_list[w],delim,vertex_type_list[j],delim,"4-chord-edge-orbit")]+=1
				elseif (v in jPath) 
					count_dict[string(vertex_type_list[v],delim,vertex_type_list[j],delim,vertex_type_list[w],delim,vertex_type_list[i],delim,"4-chord-edge-orbit")]+=1
				elseif (!(v in gamma_i) & !(v in gamma_j))
						count_dict[string(vertex_type_list[i],delim,vertex_type_list[j],delim,vertex_type_list[w],delim,vertex_type_list[v],delim,"4-tail-tri-centre-orbit")]+=1
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

		##4-path centre orbit
		for (inda,a) in enumerate(types)	
			for (indb,b) in enumerate(types[inda:end])
				if (a == b)
					#order doesn't matter here
					count_dict[string(a,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,b,delim,"4-path-centre-orbit")] = iPathTypes[inda]*jPathTypes[inda+indb-1] - count_dict[string(a,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,b,delim,"4-cycle")]
				else
					#to maintain order here, we diverge from ROssi et al and calculate each orientation separately
					count_dict[string(a,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,b,delim,"4-path-centre-orbit")] = iPathTypes[inda]*jPathTypes[inda+indb-1] - count_dict[string(a,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,b,delim,"4-cycle")]
					count_dict[string(b,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,a,delim,"4-path-centre-orbit")] = iPathTypes[inda+indb-1]*jPathTypes[inda] - count_dict[string(b,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,a,delim,"4-cycle")]
				end
	 		end
	 	end

		##4-star
		#To maintain type order here, we also have to separate. Note we enforce that the centre of the star is THIRD listed (in line with the 4-tail edge orbit layout)
		for (inda,a) in enumerate(types)	
			for (indb,b) in enumerate(types[inda:end])
				if  (a == b)
					##We only need to split here if i and j are of different types
					if (vertex_type_list[i]!=vertex_type_list[j])
						#i-centre
						count_dict[string(a,delim,b,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-star")] = 0.5*iPathTypes[inda]*(iPathTypes[inda]-1) - 	count_dict[string(a,delim,b,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-tail-edge-orbit")]

						#j-centre
						count_dict[string(a,delim,b,delim,vertex_type_list[j],delim,vertex_type_list[i],delim,"4-star")] = 0.5*jPathTypes[inda]*(jPathTypes[inda]-1) - 	count_dict[string(a,delim,b,delim,vertex_type_list[j],delim,vertex_type_list[i],delim,"4-tail-edge-orbit")]

					else # when i and j are also of same type
						count_dict[string(a,delim,b,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-star")] = 0.5*iPathTypes[inda]*(iPathTypes[inda]-1)+0.5*jPathTypes[inda]*(jPathTypes[inda]-1) - 	count_dict[string(a,delim,b,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-tail-edge-orbit")]
					end

				else # when a and b are of different types
				#again, to maintain order here, we diverge from ROssi et al and calculate each orientation separately WHEN J AND I ARE OF DIFFERENT TYPE
					if (vertex_type_list[i]!=vertex_type_list[j])
						count_dict[string(a,delim,b,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-star")] = iPathTypes[inda]*iPathTypes[inda+indb-1] - count_dict[string(a,delim,b,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-tail-edge-orbit")] - count_dict[string(b,delim,a,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-tail-edge-orbit")]#unsure if both need to be subtracted here? TEST					
						count_dict[string(a,delim,b,delim,vertex_type_list[j],delim,vertex_type_list[i],delim,"4-star")] = jPathTypes[inda]*jPathTypes[inda+indb-1] - count_dict[string(a,delim,b,delim,vertex_type_list[j],delim,vertex_type_list[i],delim,"4-tail-edge-orbit")] - count_dict[string(b,delim,a,delim,vertex_type_list[j],delim,vertex_type_list[i],delim,"4-tail-edge-orbit")]#unsure if both need to be subtracted here? TEST
					else
						
					count_dict[string(a,delim,b,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-star")] = iPathTypes[inda]*iPathTypes[inda+indb-1] +  jPathTypes[inda]*jPathTypes[inda+indb-1]  - count_dict[string(a,delim,b,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-tail-edge-orbit")] - count_dict[string(b,delim,a,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-tail-edge-orbit")]#unsure if both need to be subtracted here? TEST					
					end
				end
			end
		end	
		
		##4-tail tri-edge orbit

		for (inda,a) in enumerate(types)	
			for (indb,b) in enumerate(types[inda:end])
				if (a==b)
					##We only need to split here if i and j are of different types
					if (vertex_type_list[i]!=vertex_type_list[j])
						#i-centre
						count_dict[string(a,delim,b,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-tail-tri-edge-orbit")] = TriTypes[inda]*iPathTypes[inda] - 	count_dict[string(a,delim,b,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-chord-edge-orbit")]

						#j-centre
						count_dict[string(a,delim,b,delim,vertex_type_list[j],delim,vertex_type_list[i],delim,"4-tail-tri-edge-orbit")] = TriTypes[inda]*jPathTypes[inda] - 	count_dict[string(a,delim,b,delim,vertex_type_list[j],delim,vertex_type_list[i],delim,"4-chord-edge-orbit")]
					else # when i and j are of same type
						count_dict[string(a,delim,b,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-tail-tri-edge-orbit")] = TriTypes[inda]*iPathTypes[inda] + TriTypes[inda]*jPathTypes[inda] - count_dict[string(a,delim,b,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-chord-edge-orbit")]
					end
				else # when a and b are different types
					if (vertex_type_list[i]!=vertex_type_list[j])
						#i-centre
						count_dict[string(a,delim,b,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-tail-tri-edge-orbit")] = TriTypes[inda]*iPathTypes[inda] + TriTypes[indb]*iPathTypes[indb] - count_dict[string(a,delim,b,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-chord-edge-orbit")]- count_dict[string(b,delim,a,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-chord-edge-orbit")]

						#j-centre
						count_dict[string(a,delim,b,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-tail-tri-edge-orbit")] = TriTypes[inda]*jPathTypes[inda] + TriTypes[indb]*jPathTypes[indb] - count_dict[string(a,delim,b,delim,vertex_type_list[j],delim,vertex_type_list[i],delim,"4-chord-edge-orbit")]- count_dict[string(b,delim,a,delim,vertex_type_list[j],delim,vertex_type_list[i],delim,"4-chord-edge-orbit")]
					else # when i and j are of same type
						count_dict[string(a,delim,b,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-tri-edge-orbit")] = TriTypes[inda]*iPathTypes[inda] + TriTypes[indb]*iPathTypes[indb] + TriTypes[inda]*jPathTypes[inda] + TriTypes[indb]*jPathTypes[indb] - count_dict[string(a,delim,b,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-chord-edge-orbit")]- count_dict[string(b,delim,a,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,"4-chord-edge-orbit")]
					end
				end
			end
		end
		
		## 4-chord-centre-orbit
		for (inda,a) in enumerate(types)	
			for (indb,b) in enumerate(types[inda:end])
				if  (a == b)
					count_dict[string(a,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,b,delim,"4-chord-centre-orbit")] = 0.5*TriTypes[inda]*(TriTypes[inda]-1) - count_dict[string(a,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,b,delim,"4-clique")]
				else
					count_dict[string(a,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,b,delim,"4-chord-centre-orbit")] = TriTypes[inda]*TriTypes[indb] - count_dict[string(a,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,b,delim,"4-clique")] - count_dict[string(b,delim,vertex_type_list[i],delim,vertex_type_list[j],delim,a,delim,"4-clique")]
				end
			end
		end


#		#4-node graphlets		
#        	#find nodes that aren't connected to either i or j
#        	E = setdiff(ordered_vertices,union(gamma_j,gamma_i))
#        	##get neighbours of (neighbours of i that aren't i or neighbours of j)
#        	iiPath = setdiff.(neighbourdictfunc.(iPath),Ref(i))
#        	##get neighbours of (neighbours of j that arent j or neighbours of i)
#        	jjPath = setdiff.(neighbourdictfunc.(jPath),Ref(j))
#        	
#        	##4-PATH
#        	#4-path edge orbits, i-stem
#        	##get all neighours of i that have neighbours that are not connected to i or j
#        	istem = intersect.(iiPath,Ref(E))
#        	count_dict = add4graphlets(vertex_type_list,iPath,istem,count_dict,j,i,graphlet_type = "4-path-edge-orbit")
#        	#4-path edge orbits, j-stem
#        	#get all neighours of j that have neighbours that are not connected to i or j
#        	jstem = intersect.(jjPath,Ref(E))
#        	count_dict = add4graphlets(vertex_type_list,jPath,jstem,count_dict,i,j,graphlet_type = "4-path-edge-orbit")			
#        	#4-path centre orbits 
#        	#get all (uniquely) neighbours of j that are not connected to neighbours of i (note, this is symmetric, so do not need to do for the other way round as well)
#        	centres = Array{Array{Int,1},1}(undef,length(jPath))
#        	for p in 1:length(jPath)
#        		centres[p] = iPath[BitArray(in.(jPath[p],iiPath).*-1 .+1)]
#        	end
#        	count_dict = add4graphlets(vertex_type_list,jPath,centres,count_dict,i,j,graphlet_type = "4-path-centre-orbit")
#        
#        	## 4-STAR
#        	# 4-star istem
#        	# find pairs of neighbours of i that are not themselves connected.
#                # prioritises earlier nodes in iPath, so that later ones are not double
#                # counted.
#        	istars = Array{Array{Int,1},1}(undef,length(iPath))
#                for (ind,p) in enumerate(iPath)
#                        istars[ind] = setdiff(iPath[BitArray(in.(Ref(p),neighbourdictfunc.(iPath)).*-1 .+1)],iPath[1:ind])
#                end
#                count_dict = add4graphlets(vertex_type_list,iPath,istars,count_dict,i,j,graphlet_type = "4-star")
#        
#                #4-star jstem (as above)
#        	jstars = Array{Array{Int,1},1}(undef,length(jPath))
#                for (ind,p) in enumerate(jPath)
#                        jstars[ind] = setdiff(jPath[BitArray(in.(Ref(p),neighbourdictfunc.(jPath)).*-1 .+1)],jPath[1:ind])
#                end
#                count_dict = add4graphlets(vertex_type_list,jPath,jstars,count_dict,j,i,graphlet_type = "4-star")
#        
#        	## 4-TAIL
#        	#4-tail-edge-orbit istem
#        	#Similar to 4-star, but instead want those neighbour pairs that ARE connected
#        	itails = Array{Array{Int,1},1}(undef,length(iPath))
#                for (ind,p) in enumerate(iPath)
#                        itails[ind] = setdiff(iPath[in.(Ref(p),neighbourdictfunc.(iPath))],iPath[1:ind])
#                end
#                count_dict = add4graphlets(vertex_type_list,iPath,itails,count_dict,i,j,graphlet_type = "4-tail-edge-orbit")
#        	
#		
#
#                #4-tail-edge-orbit jstem (as above)
#        	jtails = Array{Array{Int,1},1}(undef,length(jPath))
#                for (ind,p) in enumerate(jPath)
#                        jtails[ind] = setdiff(jPath[in.(Ref(p),neighbourdictfunc.(jPath))],jPath[1:ind])
#                end
#                count_dict = add4graphlets(vertex_type_list,jPath,jtails,count_dict,j,i,graphlet_type = "4-tail-edge-orbit")
#        	
#		#4-tail-tri-edge
#		##this is instead based off protrusions off the triangles of i and j
#		#first we look at protrusions off of i:
#        	ittails = Array{Array{Int,1},1}(undef,length(iPath))
#		##find neighbours of i that are not neighbours of j (iPath) that are also not neighbours of k (Tri)
#		for (ind,p) in enumerate(iPath)
#                        ittails[ind] = Tri[BitArray(in.(Ref(p),neighbourdictfunc.(Tri)).*-1 .+1)]
#                end
#                count_dict = add4graphlets(vertex_type_list,iPath,ittails,count_dict,i,j,graphlet_type = "4-tail-tri-edge-orbit")
#		
#		#the same for j:
#        	jttails = Array{Array{Int,1},1}(undef,length(jPath))
#		##find neighbours of j that are not neighbours of i (jPath) that are also not neighbours of k (Tri)
#		for (ind,p) in enumerate(jPath)
#                        jttails[ind] = Tri[BitArray(in.(Ref(p),neighbourdictfunc.(Tri)).*-1 .+1)]
#                end
#                count_dict = add4graphlets(vertex_type_list,jPath,jttails,count_dict,j,i,graphlet_type = "4-tail-tri-edge-orbit")
#
#		#4-tail-tri-centre
#		## find neighbours of the third node in the triangle that are not connected to i or j
#                kttails = intersect.(setdiff.(neighbourdictfunc.(Tri),Ref([i,j])),Ref(E))		
#                count_dict = add4graphlets(vertex_type_list,Tri,kttails,count_dict,i,j,graphlet_type = "4-tail-tri-centre-orbit")
#		
#
#
#        	## 4-CYCLE
#                #
#                #Again, this is symmetric, so we only need to count this from one perspective
#                cycles = Array{Array{Int,1},1}(undef,length(iPath))
#                for (ind,p) in enumerate(iPath)
#                        cycles[ind] = jPath[in.(Ref(p),neighbourdictfunc.(jPath))]
#                end
#                count_dict = add4graphlets(vertex_type_list,iPath,cycles,count_dict,i,j,graphlet_type = "4-cycle")
#                
#  		##4-CHORD
#		#
#		#4-chord, i edge
#		#similar to 4-tail-tri-centre, but we instead find neighbours of k that are connected to one of i or j
#		## find neighbours of the third node in the triangle that are not connected to i or j
#		ichord = intersect.(setdiff.(neighbourdictfunc.(Tri),Ref([i,j])),Ref(iPath))		
#                count_dict = add4graphlets(vertex_type_list,Tri,ichord,count_dict,i,j,graphlet_type = "4-chord-edge-orbit")
#		
#		jchord = intersect.(setdiff.(neighbourdictfunc.(Tri),Ref([i,j])),Ref(jPath))		
#                count_dict = add4graphlets(vertex_type_list,Tri,jchord,count_dict,j,i,graphlet_type = "4-chord-edge-orbit")
#
#		##4-chord centre orbit
#		#now we need to identify triangles (k in Tri) that don't share an edge with one another
#        	chord = Array{Array{Int,1},1}(undef,length(Tri))
#                for (ind,p) in enumerate(Tri)
#                        chord[ind] = setdiff(Tri[BitArray(in.(Ref(p),neighbourdictfunc.(Tri)).*-1 .+1)],Tri[1:ind])
#                end
#                count_dict = add4graphlets(vertex_type_list,Tri,chord,count_dict,i,j,graphlet_type = "4-chord-centre-orbit")
#		
#		##4-CLIQUE
#		#Finally, we need to identify triangles (k in Tri) that DO share an edge with one another
#        	clique = Array{Array{Int,1},1}(undef,length(Tri))
#                 for (ind,p) in enumerate(Tri)
#                        clique[ind] = setdiff(Tri[in.(Ref(p),neighbourdictfunc.(Tri))],Tri[1:ind])
#                end
#                count_dict = add4graphlets(vertex_type_list,Tri,clique,count_dict,i,j,graphlet_type = "4-clique")
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
	for h in 1 :size(edgelist,1)
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
				#we only want to switch if interior needs switching
				if (graphlet_names[el][[2,3]][1]!=sort(graphlet_names[el][[2,3]])[1])
					graphlet_names[el][[2,3]] = sort(graphlet_names[el][[2,3]])
					graphlet_names[el][[1,4]] = graphlet_names[el][[4,1]]
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
