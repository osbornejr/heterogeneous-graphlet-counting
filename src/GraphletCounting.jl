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
	list=fill("",length(nodelist))
	delim="_"
	for (ind ,n) in enumerate(nodelist)
		list[ind]=string(vertexlist[i],delim,vertexlist[j],delim,vertexlist[n],delim,graphlet_type)
	end
	#list=vertexlist[i].*"_".*vertexlist[j].*"_".*vertexlist[nodelist].*"_".*graphlet_type	
	x=Dict{String,Integer}([(i,count(x->x==i,list)) for i in unique(list)]) 
	merge!(+,count_dict,x)
	return count_dict
end

#### RT-POL specific archive
##import data
#import CSV
#edgelist=CSV.read("data/rt-pol/soc-political-retweet.edges",header=["V1","V2","ref"],type=Int)
#edgelist=edgelist[:,1:2]
#vertexlist=CSV.read("data/rt-pol/soc-political-retweet.node_labels",header=["V","type"],type=Int)

#fix vertex values to iterate from 1, not 0
#edgelist[:,1:2]=edgelist[:,1:2].+1
#vertexlist.V=vertexlist.V.+1

#change label names
#vertexlist.type=replace(vertexlist.type,Pair(1,"right"),Pair(2,"left"))

#create graph from edgelist (only really needed for visualisation, moving on now)
#using LightGraphs
#using LightGraphs.SimpleGraphs
#using GraphPlot
#g₁=SimpleGraph(map(SimpleEdge,Tuple.(eachrow(edgelist))))
#gplot(g₁) 
####


#start writing required functions ad hoc
#change to arrays for performance (also remove multi-edges from edgelist here for now)
#edgelist=unique(sort(convert(Array{Int,2},edgelist),dims=2),dims=1)

function count_graphlets(vertex_type_list::Array{String,1},edgelist::Array{Pair,1})
	
	#vertexlist=convert(Array{String},vertexlist[:,2]);
	vertexlist = vertex_type_list
	#get neighbourhood for each vertex in advance (rather than calling per-edge)
	neighbourdict=Neighbours(edgelist)
	#preallocate array to store each edge's graphlet dictionary 
	Chi=Array{Dict}(undef,size(edgelist,1));
	
	#per edge process
	for h in 1:size(edgelist,1)
		count_dict = Dict{String,Int}()
	#	h=1	
		i = edgelist[h].first
		j = edgelist[h].second
		gamma_i = neighbourdict[i]	
		gamma_j = neighbourdict[j]
		
		#three node graphlets
		Tri = intersect(gamma_i,gamma_j)
		count_dict = add3graphlets(vertexlist,Tri,count_dict,i,j,graphlet_type="3-tri")
		jPath = setdiff(gamma_j,union(gamma_i,i))
		count_dict = add3graphlets(vertexlist,jPath,count_dict,i,j,graphlet_type="3-path")
		iPath = setdiff(gamma_i,union(gamma_j,j))
		count_dict = add3graphlets(vertexlist,iPath,count_dict,j,i,graphlet_type="3-path")
		
		#Save count dictionary for this edge
		Chi[h] = count_dict
	end
	#total counts for each graphlet
	total_counts = reduce(mergecum,Chi)
	
	#reorder names to merge orbits
	graphlet_names = (split.(collect(keys(total_counts)),"_"))
	for el in 1:size(graphlet_names,1)
		#do not reorder for x-y-x paths (different orbit to other 3-paths) 
		if !(graphlet_names[el][1]!=graphlet_names[el][2] && graphlet_names[el][1]==graphlet_names[el][3] && graphlet_names[el][4]=="3-path")
	          graphlet_names[el][1:3]=sort(graphlet_names[el][1:3])
	        end
	end
	orbit_counts = Dict{String,Int}()
	orbits = join.(graphlet_names,"_")
	count_values = collect(values(total_counts))
	for orb in 1:size(orbits,1)
		orbit_counts[orbits[orb]] = get(orbit_counts,orbits[orb],0)+count_values[orb]	
	end
	return orbit_counts
end


#count triangles for comparison
#tri_counts=[orbit_counts[i] for i in collect(keys(orbit_counts)) if occursin("3-tri",i)]
#tri_counts=tri_counts./sum(tri_counts)

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
	graphlet_counts = count_graphlets(vertexlist[:,2],edgelist)
	
	#Statistical significance
	zscores = Dict{String,Float64}()
	for g in collect(keys(graphlet_counts))
		rand_vals = filter(:graphlet=>x->x==g,null_model_df)[!,:value]
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
		draw(SVG(plotfile),p)
		print("Stat plot saved to $plotfile.")
	end
	#null_model_sum = reduce(mergecum,null_model_calc)
	return zscores
end	
