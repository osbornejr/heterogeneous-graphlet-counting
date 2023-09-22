
"""

    graphlet_counts_per_community

Find graphlet counts for each community subgraph in network.

"""
function graphlet_counts_per_community(vertexlist::Vector{<:AbstractString},edgelist::Union{Array{Pair{Int,Int},1},Array{Pair,1}},community_groups::Vector{Int},graphlet_size::Int=3;run_method::String="serial",progress::Bool=false,recursive::Bool=true)
    
    if(length(vertexlist)!==length(community_groups))
        throw(ArgumentError("vector of community groups must be same length as vertex type list"))
    end
    
    comms = unique(community_groups)
    ##store counts for each comm here
    Comm_arr = Dict{String,Int}[]
    for c in comms 
        ##get sub vertexlist (using index as ref)
        c_v = findall(.==(c),community_groups)
        ##get sub edgelist
        c_e = edgelist[findall(.==(2),sum.(is_community_edge.(edgelist,Ref(c_v))))]
        c_g = GraphletCounting.count_graphlets(vertexlist,c_e,graphlet_size,run_method=run_method,progress=progress)           
        push!(Comm_arr,c_g)
    end
    return Comm_arr
end

"""
    
    is_community_edge

Determine whether the vertices of an edge pair are both within, partially within or without the community vertexlist. 

"""
function is_community_edge(edge::Pair,vertexlist::Vector{Int})
    return (in(first(edge),vertexlist),in(last(edge),vertexlist))
end

"""
    
    cross_community_edges

Find all edges that exist between communities for a given commmunity partition.  

"""
function cross_community_edges(edgelist::Vector{Pair},community_partition::Vector{Int})
    
    ##get all unique communities
    comms = unique(community_partition)
    
    ##variables
    n_comms = length(comms)
    n_edges = length(edgelist)

    ##store edge status for each community
    edge_arr = Vector{Int}[]
    for c in comms 
        ##get sub vertexlist (using index as ref)
        c_v = findall(.==(c),community_partition)
        ##get edge status for this community
        c_e = sum.(is_community_edge.(edgelist,Ref(c_v)))
        #store
        push!(edge_arr,c_e)
    end
    ## check status of each edge over all communities
    comb = hcat(edge_arr...)
    edge_cat = [max(comb[x,:]...) for x in 1:length(edgelist)]
    edge_bool = (edge_cat.-2).*-1
    return BitVector(edge_bool)
    ## get more detailed per community edge details      
    # first isolate in and cross community edge matrices  
    #in_comm =(comb.-1).>0
    #cross_comm = (comb.-1).==0
    ##dummy matrix with repeated vector of comm ids 
    #comm_matrix = reshape(repeat(collect(1:n_comms),n_edges),n_comms,:)'
    ## vector corresponding to community of in community edges 
    #ins = sum(comm_matrix.*in_comm,dims=2) 
    # next find communities associated with cross edges 
    #comm_matrix.*cross_comm
    #return ins
end

##make dictionaries callable/broadcastable
(d::Dict)(k) = d[k]

"""
    edgelist_community_status

Give status of each edge in terms of the given community partition.
"""
function edgelist_community_status(edgelist::Vector{Pair},community_partition::Vector{Int})
   # dict to match vertex to its community
   comm_dict = Dict(Pair.(1:length(community_partition),community_partition))
   ##call dict to get community edgelist 
   comm_edgelist = Pair.(comm_dict.(first.(edgelist)),comm_dict.(last.(edgelist)))
   return comm_edgelist
end


"""
    types_per_community

Find type distribution for each community subgraph in network.

"""
function types_per_community(vertexlist::Vector{<:AbstractString},community_groups::Vector{Int})
    types = unique(vertexlist)
    ##first print out percentage for whole network
    for t in types
        print("In whole network: $(round(sum(vertexlist.==t)/length(vertexlist),digits=4)) $(t)\n")
    end
    print("\n")
    ## then for each community
    comms = unique(community_groups)
    for c in comms
        c_vl = vertexlist[community_groups.==c] 
        for t in types
            print("In community $(c): $(round(sum(c_vl.==t)/length(c_vl),digits=4)) $(t)\n")

        end
    print("\n")
    end
end


