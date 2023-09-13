
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
Determine whether the vertices of an edge pair are both within, partially within or without the community. 

"""
function is_community_edge(edge::Pair,vertexlist::Vector{Int})
    return (in(first(edge),vertexlist),in(last(edge),vertexlist))
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


