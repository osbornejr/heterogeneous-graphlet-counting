
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

#function add4graphlets(typelist::Array{String,1},nodelist1::Array{Int,1},nodelist2::Array{Array{Int,1},1},count_dict::Dict{String,Int},i::Int,j::Int;graphlet_type::String,order::Array{Int,1}=[1,2,3,4])
#
#	#check to see if candidate lists are empty
#	if length(nodelist1)==0 
#		return count_dict
#	elseif sum(length.(nodelist2))==0
#		return count_dict
#	else ##we do have some graphlets to add!
#		delim = "_"
#		for (ind,n) in enumerate(nodelist1)
#			#check if specific list is empty
#			if length(nodelist2[ind])==0
#				x = Dict{String,Int}()
#			else
#				list = fill("",sum(length(nodelist2[ind])))
#				if (graphlet_type=="4-path-edge-orbit")
#					for (indd,m) in enumerate(nodelist2[ind])
#						#ordering i_j_jneigh_j_neighneigh
#						list[indd] = string(typelist[i],delim,typelist[j],delim,typelist[n],delim,typelist[m],delim,graphlet_type)
#					end
#				end
#				if (graphlet_type == "4-path-centre-orbit")
#					for (indd,m) in enumerate(nodelist2[ind])
#						#ordering neigh_i_i_j_jneigh
#						list[indd] = string(typelist[m],delim,typelist[i],delim,typelist[j],delim,typelist[n],delim,graphlet_type)
#					end
#				end
#                                if (graphlet_type == "4-tail-edge-orbit" )
#					for (indd,m) in enumerate(nodelist2[ind])
#						#ordering neighi_neighi_i_j
#                                                list[indd] = string(typelist[n],delim,typelist[m],delim,typelist[i],delim,typelist[j],delim,graphlet_type)
#					end
#				end
#                                if (graphlet_type == "4-star")
#					for (indd,m) in enumerate(nodelist2[ind])
#						#ordering j_neighi_neighi_i
#                                                list[indd] = string(typelist[j],delim,typelist[n],delim,typelist[m],delim,typelist[i],delim,graphlet_type)
#					end
#                                        
#                                end
#                                if (graphlet_type == "4-tail-tri-edge-orbit" )
#					for (indd,m) in enumerate(nodelist2[ind])
#						#ordering k_j_i_neigh_i
#						list[indd] = string(typelist[m],delim,typelist[j],delim,typelist[i],delim,typelist[n],delim,graphlet_type)
#					end
#                                        
#                                end
#                                if (graphlet_type == "4-tail-tri-centre-orbit" )
#					for (indd,m) in enumerate(nodelist2[ind])
#						#ordering i_j_k_neigh_k
#						list[indd] = string(typelist[i],delim,typelist[j],delim,typelist[n],delim,typelist[m],delim,graphlet_type)
#					end
#                                        
#                                end
#                                if (graphlet_type == "4-cycle")
#					for (indd,m) in enumerate(nodelist2[ind])
#						#ordering neighi_i_j_neighj
#                                                list[indd] = string(typelist[n],delim,typelist[i],delim,typelist[j],delim,typelist[m],delim,graphlet_type)
#					end
#                                        
#                                end
#                                if (graphlet_type == "4-chord-edge-orbit")
#					for (indd,m) in enumerate(nodelist2[ind])
#						#ordering neighik_i_k_j
#                                                list[indd] = string(typelist[m],delim,typelist[i],delim,typelist[n],delim,typelist[j],delim,graphlet_type)
#					end
#                                        
#                                end
#                                if (graphlet_type == "4-chord-centre-orbit")
#					for (indd,m) in enumerate(nodelist2[ind])
#						#ordering neighij_i_j_k
#                                                list[indd] = string(typelist[m],delim,typelist[i],delim,typelist[j],delim,typelist[n],delim,graphlet_type)
#					end
#                                        
#                                end
#                                if (graphlet_type == "4-clique")
#					for (indd,m) in enumerate(nodelist2[ind])
#						#ordering k_i_j_neigh_kij
#                                                list[indd] = string(typelist[n],delim,typelist[i],delim,typelist[j],delim,typelist[m],delim,graphlet_type)
#					end
#                                        
#                                end
#				x = Dict{String,Int}([(p,count(x->x==p,list)) for p in unique(list)])
#			end
#		merge!(+,count_dict,x)
#		end
#	end
#	return count_dict
#end

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
						#ordering neighi_neighi_i_j
						count_dict[string(typelist[n],delim,typelist[m],delim,typelist[i],delim,typelist[j],delim,graphlet_type)]+=1
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

function per_edge_counts_old(edge::Int,vertex_type_list::Array{String,1},edgelist::Union{Array{Pair{Int,Int},1},Array{Pair,1}},graphlet_size::Int,neighbourdict::Dict{Int,Vector{Int}},neighbourdictfunc::Function,ordered_vertices::Array{Int,1})
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
		#4-node graphlets		
        	#find nodes that aren't connected to either i or j
        	E = setdiff(ordered_vertices,union(gamma_j,gamma_i))
        	##get neighbours of (neighbours of i that aren't i or neighbours of j)
        	iiPath = setdiff.(neighbourdictfunc.(iPath),Ref(i))
        	##get neighbours of (neighbours of j that arent j or neighbours of i)
        	jjPath = setdiff.(neighbourdictfunc.(jPath),Ref(j))
        	
        	##4-PATH
        	#4-path edge orbits, i-stem
        	##get all neighours of i that have neighbours that are not connected to i or j
        	istem = intersect.(iiPath,Ref(E))
        	count_dict = add4graphlets(vertex_type_list,iPath,istem,count_dict,j,i,graphlet_type = "4-path-edge-orbit")
        	#4-path edge orbits, j-stem
        	#get all neighours of j that have neighbours that are not connected to i or j
        	jstem = intersect.(jjPath,Ref(E))
        	count_dict = add4graphlets(vertex_type_list,jPath,jstem,count_dict,i,j,graphlet_type = "4-path-edge-orbit")			
        	#4-path centre orbits 
        	#get all (uniquely) neighbours of j that are not connected to neighbours of i (note, this is symmetric, so do not need to do for the other way round as well)
        	centres = Array{Array{Int,1},1}(undef,length(jPath))
        	for p in 1:length(jPath)
        		centres[p] = iPath[BitArray(in.(jPath[p],iiPath).*-1 .+1)]
        	end
        	count_dict = add4graphlets(vertex_type_list,jPath,centres,count_dict,i,j,graphlet_type = "4-path-centre-orbit")
        
        	## 4-STAR
        	# 4-star istem
        	# find pairs of neighbours of i that are not themselves connected.
                # prioritises earlier nodes in iPath, so that later ones are not double
                # counted.
        	istars = Array{Array{Int,1},1}(undef,length(iPath))
                for (ind,p) in enumerate(iPath)
                        istars[ind] = setdiff(iPath[BitArray(in.(Ref(p),neighbourdictfunc.(iPath)).*-1 .+1)],iPath[1:ind])
                end
                count_dict = add4graphlets(vertex_type_list,iPath,istars,count_dict,i,j,graphlet_type = "4-star")
        
                #4-star jstem (as above)
        	jstars = Array{Array{Int,1},1}(undef,length(jPath))
                for (ind,p) in enumerate(jPath)
                        jstars[ind] = setdiff(jPath[BitArray(in.(Ref(p),neighbourdictfunc.(jPath)).*-1 .+1)],jPath[1:ind])
                end
                count_dict = add4graphlets(vertex_type_list,jPath,jstars,count_dict,j,i,graphlet_type = "4-star")
        
        	## 4-TAIL
        	#4-tail-edge-orbit istem
        	#Similar to 4-star, but instead want those neighbour pairs that ARE connected
        	itails = Array{Array{Int,1},1}(undef,length(iPath))
                for (ind,p) in enumerate(iPath)
                        itails[ind] = setdiff(iPath[in.(Ref(p),neighbourdictfunc.(iPath))],iPath[1:ind])
                end
                count_dict = add4graphlets(vertex_type_list,iPath,itails,count_dict,i,j,graphlet_type = "4-tail-edge-orbit")
        	
		

                #4-tail-edge-orbit jstem (as above)
        	jtails = Array{Array{Int,1},1}(undef,length(jPath))
                for (ind,p) in enumerate(jPath)
                        jtails[ind] = setdiff(jPath[in.(Ref(p),neighbourdictfunc.(jPath))],jPath[1:ind])
                end
                count_dict = add4graphlets(vertex_type_list,jPath,jtails,count_dict,j,i,graphlet_type = "4-tail-edge-orbit")
        	
		#4-tail-tri-edge
		##this is instead based off protrusions off the triangles of i and j
		#first we look at protrusions off of i:
        	ittails = Array{Array{Int,1},1}(undef,length(iPath))
		##find neighbours of i that are not neighbours of j (iPath) that are also not neighbours of k (Tri)
		for (ind,p) in enumerate(iPath)
                        ittails[ind] = Tri[BitArray(in.(Ref(p),neighbourdictfunc.(Tri)).*-1 .+1)]
                end
                count_dict = add4graphlets(vertex_type_list,iPath,ittails,count_dict,i,j,graphlet_type = "4-tail-tri-edge-orbit")
		
		#the same for j:
        	jttails = Array{Array{Int,1},1}(undef,length(jPath))
		##find neighbours of j that are not neighbours of i (jPath) that are also not neighbours of k (Tri)
		for (ind,p) in enumerate(jPath)
                        jttails[ind] = Tri[BitArray(in.(Ref(p),neighbourdictfunc.(Tri)).*-1 .+1)]
                end
                count_dict = add4graphlets(vertex_type_list,jPath,jttails,count_dict,j,i,graphlet_type = "4-tail-tri-edge-orbit")

		#4-tail-tri-centre
		## find neighbours of the third node in the triangle that are not connected to i or j
                kttails = intersect.(setdiff.(neighbourdictfunc.(Tri),Ref([i,j])),Ref(E))		
                count_dict = add4graphlets(vertex_type_list,Tri,kttails,count_dict,i,j,graphlet_type = "4-tail-tri-centre-orbit")
		


        	## 4-CYCLE
                #
                #Again, this is symmetric, so we only need to count this from one perspective
                cycles = Array{Array{Int,1},1}(undef,length(iPath))
                for (ind,p) in enumerate(iPath)
                        cycles[ind] = jPath[in.(Ref(p),neighbourdictfunc.(jPath))]
                end
                count_dict = add4graphlets(vertex_type_list,iPath,cycles,count_dict,i,j,graphlet_type = "4-cycle")
                
  		##4-CHORD
		#
		#4-chord, i edge
		#similar to 4-tail-tri-centre, but we instead find neighbours of k that are connected to one of i or j
		## find neighbours of the third node in the triangle that are not connected to i or j
		ichord = intersect.(setdiff.(neighbourdictfunc.(Tri),Ref([i,j])),Ref(iPath))		
                count_dict = add4graphlets(vertex_type_list,Tri,ichord,count_dict,i,j,graphlet_type = "4-chord-edge-orbit")
		
		jchord = intersect.(setdiff.(neighbourdictfunc.(Tri),Ref([i,j])),Ref(jPath))		
                count_dict = add4graphlets(vertex_type_list,Tri,jchord,count_dict,j,i,graphlet_type = "4-chord-edge-orbit")

		##4-chord centre orbit
		#now we need to identify triangles (k in Tri) that don't share an edge with one another
        	chord = Array{Array{Int,1},1}(undef,length(Tri))
                for (ind,p) in enumerate(Tri)
                        chord[ind] = setdiff(Tri[BitArray(in.(Ref(p),neighbourdictfunc.(Tri)).*-1 .+1)],Tri[1:ind])
                end
                count_dict = add4graphlets(vertex_type_list,Tri,chord,count_dict,i,j,graphlet_type = "4-chord-centre-orbit")
		
		##4-CLIQUE
		#Finally, we need to identify triangles (k in Tri) that DO share an edge with one another
        	clique = Array{Array{Int,1},1}(undef,length(Tri))
                for (ind,p) in enumerate(Tri)
                        clique[ind] = setdiff(Tri[in.(Ref(p),neighbourdictfunc.(Tri))],Tri[1:ind])
                end
                count_dict = add4graphlets(vertex_type_list,Tri,clique,count_dict,i,j,graphlet_type = "4-clique")
 	end
	return count_dict 
end

function count_graphlets_old(vertex_type_list::Array{String,1},edgelist::Union{Array{Pair{Int,Int},1},Array{Pair,1}},graphlet_size::Int=3)


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
	for h in 1:size(edgelist,1)
		count_dict = per_edge_counts_old(h,vertex_type_list,edgelist,graphlet_size,neighbourdict,neighbourdictfunc,ordered_vertices)        
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
