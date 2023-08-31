
"""

    graphlet_counts_per_community

Find graphlet counts for each community subgraph in network.

"""
function graphlet_counts_per_community(vertexlist::Vector{<:AbstractString},edgelist::Union{Array{Pair{Int,Int},1},Array{Pair,1}},community_groups::Vector{Int},graphlet_size::Int=3;run_method::String="serial",progress::Bool=false,recursive::Bool=true)
  

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


