using DataStructures,Distributed,ProgressMeter,StatsBase

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

function graphlet_string(a::AbstractString,b::AbstractString,c::AbstractString,graphlet::AbstractString,delim::AbstractString)
    return string(a,delim,b,delim,c,delim,graphlet)
end
function graphlet_string(a::AbstractString,b::AbstractString,c::AbstractString,d::AbstractString,graphlet::AbstractString,delim::AbstractString)
    return string(a,delim,b,delim,c,delim,d,delim,graphlet)
end


function per_edge_relationships(edge::Int,vertex_type_list::Vector{<:AbstractString},edgelist::Union{Array{Pair{Int,Int},1},Array{Pair,1}},graphlet_size::Int,neighbourdict::Dict{Int,Vector{Int}},outfile::String;write_out=false)
    h=edge  
#   # get nodes i and j for this edge
        i = edgelist[h].first
        j = edgelist[h].second
        #get neighbourhoods of i and j
        gamma_i = neighbourdict[i]  
        gamma_j = neighbourdict[j]
        
        delim = "_"
        #three node graphlets
    ##more efficient loop based method
    #rel = DefaultDict{Int,Int}(0)
    rel = zeros(Int,length(vertex_type_list))
    Tri = Array{Int,1}()
    iPath = Array{Int,1}()
    jPath = Array{Int,1}()
    #prelim cycle through i neighbours
    for k in gamma_i
        if (k!=j)
            rel[k] = 1
        end
    end
    for k in gamma_j
        if(k!=i)
            if (rel[k]==1)
                ##triangle
                append!(Tri,k)
                rel[k] = 3
            else
                #j-path
                append!(jPath,k)
                rel[k] = 2
            end
        end
    end
        
    for k in gamma_i
        if (rel[k]==1)
            #ipaths
            append!(iPath,k)
        end
    end 
        if (graphlet_size==4)
            #matrix to store 4 node relationships in; each row corresponds to the relevant entry in rel TODO might be better for memory to only have rows associated with
            #nonzero entries in rel?
            Rel = zeros(Int,length(vertex_type_list),length(vertex_type_list))
            ##for reference, relationships are coded as follows:
            #1 = 3-path i centre
            #2 = 3-path j centre
            #3 = 3-tri
            #4 = 4-path
            #5 = 4-star
            #6 = 4-tail
            #7 = 4-cycle
            #8 = 4-chord iedge
            #9 = 4-chord jedge
            #10 = 4-clique
            #The specific orbits of each four node graphlet should be recoverable from the first column in Rel (1,2 or 3) but otherwise we might need more relationship types 
            #in the fournode case
            for w in  iPath
                for v in neighbourdict[w]
                    if (v==i)

                    elseif (!(v in gamma_i) & !(v in gamma_j))
                        Rel[w,v] = 4
                    elseif ((v in iPath) & (v < w))
                        Rel[w,v] = 6
                    end

                end
            end
                        
            for w in jPath 
                for v in neighbourdict[w]
                    if (v==j)
                        #do nothing
                    elseif (!(v in gamma_i) & !(v in gamma_j))
                        Rel[w,v] = 4
                    elseif ((v in jPath) & (v < w))
                        Rel[w,v] = 6
                    elseif (v in iPath)
                        Rel[w,v] = 7
                    end

                end
            end
    
            for w in Tri 
                for v  in neighbourdict[w]
                    if (v==i|v==j)
                        #do nothing
                    elseif ((v in Tri) & (v < w)) 
                        Rel[w,v] = 10
                        ## separating the processes here so that we can maintain the right type ordering 
                    elseif (v in iPath) 
                        Rel[w,v] = 8
                    elseif (v in jPath) 
                        Rel[w,v] = 9
                    elseif (!(v in gamma_i) & !(v in gamma_j))
                        Rel[w,v] = 6
                    end             
                end
            end
        
        end 
    
    ###extract relationships vector from matrix... I think it is better to do here as storing matrix for each edge will surely eat up all memory.
    begin
        ships = Array{Tuple{Int,Int,Int,Int,String},1}()
        #ipaths
        append!(ships,[(0,j,i,x,"3-path") for x in findall(==(1),rel)])
        #jpaths
        append!(ships,[(0,i,j,x,"3-path") for x in findall(==(2),rel)])
        #triangles
        append!(ships,[(0,i,j,x,"3-tri") for x in findall(==(3),rel)])
        
        if (graphlet_size==4)
            #4-paths iedge
            append!(ships,[(j,i,x,y,"4-path") for x in findall(==(1),rel) for y in findall(==(4),Rel[x,:])])
            #4-paths jedge
            append!(ships,[(i,j,x,y,"4-path") for x in findall(==(2),rel) for y in findall(==(4),Rel[x,:])])


            #4-tails icentre
            append!(ships,[(y,x,i,j,"4-tail") for x in findall(==(1),rel) for y in findall(==(6),Rel[x,:])])
            #4-tails jcentre
            append!(ships,[(y,x,j,i,"4-tail") for x in findall(==(2),rel) for y in findall(==(6),Rel[x,:])])
            #4-tails tricentre
            append!(ships,[(i,j,x,y,"4-tail") for x in findall(==(3),rel) for y in findall(==(6),Rel[x,:])])

            #4-cycles
            append!(ships,[(i,j,x,y,"4-cycle") for x in findall(==(2),rel) for y in findall(==(7),Rel[x,:])])

            #4-chord iedge orbit
            append!(ships,[(j,i,x,y,"4-chord") for x in findall(==(3),rel) for y in findall(==(8),Rel[x,:])])
            #4-chord jedge orbit
            append!(ships,[(i,j,x,y,"4-chord") for x in findall(==(3),rel) for y in findall(==(9),Rel[x,:])])

            #4-clique
            append!(ships,[(x,i,j,y,"4-clique") for x in findall(==(3),rel) for y in findall(==(10),Rel[x,:])])

            ##combinatorials might be more difficult...
            #4-paths centre orbit
            append!(ships,[(x,i,j,y,"4-path") for x in findall(==(1),rel) for y in findall(==(2),rel) if Rel[y,x]!=7])

            #4-chords centre orbit
            append!(ships,[(x,i,j,y,"4-chord") for x in findall(==(3),rel) for y in findall(==(3),rel) if (y<x && Rel[x,y]!=10)])

            #4-stars i centre
            append!(ships,[(x,j,i,y,"4-star") for x in findall(==(1),rel) for y in findall(==(1),rel) if (y<x && Rel[x,y]!=6)])

            #4-stars j centre
            append!(ships,[(x,i,j,y,"4-star") for x in findall(==(2),rel) for y in findall(==(2),rel) if (y<x && Rel[x,y]!=6)])

            #4-tails tri i edge
            append!(ships,[(x,j,i,y,"4-tail") for x in findall(==(3),rel) for y in findall(==(1),rel) if (Rel[x,y]!=8)])

            #4-tails tri j edge
            append!(ships,[(x,i,j,y,"4-tail") for x in findall(==(3),rel) for y in findall(==(2),rel) if (Rel[x,y]!=9)])
        end
    end
    
    if (write_out)
        ##write to file
        file = "$outfile/relationships_$(myid()).csv"
        open(file,"a",lock=true) do io
            for line in ships
                write(io,"$(join([line...],","))\n")
            end
            #ticket =false
        end
        ##return ticket to channel
        return nothing
    else
        return ships
    end
end

function per_edge_relationships_alt(edge::Int,vertex_type_list::Vector{<:AbstractString},edgelist::Union{Array{Pair{Int,Int},1},Array{Pair,1}},graphlet_size::Int,neighbourdict::Dict{Int,Vector{Int}})
    ###OUTDATED
    throw(MethodError("method is depreceated"))
    ###Turns out, this method is slower than old method, with seemingly no memory gains! strange, but it might prove better on LARGER networks (if they are ever possible memory wise).
    ### fix to clear "sticky" memory
    GC.gc()
    ccall(:malloc_trim, Cvoid, (Cint,), 0)
    h=edge  
#   # get nodes i and j for this edge
        i = edgelist[h].first
        j = edgelist[h].second
        #get neighbourhoods of i and j
        gamma_i = neighbourdict[i]  
        gamma_j = neighbourdict[j]
        
        delim = "_"
        #three node graphlets
    ##more efficient loop based method
    #rel = DefaultDict{Int,Int}(0)
    rel = zeros(Int,length(vertex_type_list))
    Tri = Array{Int,1}()
    iPath = Array{Int,1}()
    jPath = Array{Int,1}()
    #prelim cycle through i neighbours
    for k in gamma_i
        if (k!=j)
            rel[k] = 1
        end
    end
    for k in gamma_j
        if(k!=i)
            if (rel[k]==1)
                ##triangle
                append!(Tri,k)
                rel[k] = 3
            else
                #j-path
                append!(jPath,k)
                rel[k] = 2
            end
        end
    end
        
    for k in gamma_i
        if (rel[k]==1)
            #ipaths
            append!(iPath,k)
        end
    end 
        if (graphlet_size==4)
        #Rel = zeros(Int,length(vertex_type_list),length(vertex_type_list))
        #For efficiency, store relationship info in dict here rather than large matrix
        Rel =  DefaultDict{Int,Vector{Pair}}(Vector{Pair{Int,Int}}())
        ##for reference, relationships are coded as follows:
        #1 = 3-path i centre
        #2 = 3-path j centre
        #3 = 3-tri
        #4 = 4-path
        #5 = 4-star
        #6 = 4-tail
        #7 = 4-cycle
        #8 = 4-chord iedge
        #9 = 4-chord jedge
        #10 = 4-clique
        #The specific orbits of each four node graphlet should be recoverable from the first column in Rel (1,2 or 3) but otherwise we might need more relationship types 
        #in the fournode case
        for w in  iPath
            for v in neighbourdict[w]
                if (v==i)

                elseif (!(v in gamma_i) & !(v in gamma_j))
                        push!(Rel[4],w=>v)
                elseif ((v in iPath) & (v < w))
                        push!(Rel[6],w=>v)
                end

            end
        end
                        
        for w in jPath 

            for v in neighbourdict[w]
                if (v==j)
                #do nothing
                elseif (!(v in gamma_i) & !(v in gamma_j))
                        push!(Rel[4],w=>v)
                elseif ((v in jPath) & (v < w))
                        push!(Rel[6],w=>v)
                elseif (v in iPath)
                        push!(Rel[7],w=>v)
                end

            end
        end
    
        for w in Tri 
            for v  in neighbourdict[w]
                if (v==i|v==j)
                #do nothing
                elseif ((v in Tri) & (v < w)) 
                        push!(Rel[10],w=>v)
                ## separating the processes here so that we can maintain the right type ordering 
                elseif (v in iPath) 
                        push!(Rel[8],w=>v)
                elseif (v in jPath) 
                        push!(Rel[9],w=>v)
                elseif (!(v in gamma_i) & !(v in gamma_j))
                        push!(Rel[6],w=>v)
                end             
            end
        end
    
    ###extract relationships vector from matrix... I think it is better to do here as storing matrix for each edge will surely eat up all memory.
    ships = Array{Tuple{Int,Int,Int,Int,String},1}()
    #ipaths
    append!(ships,[(0,j,i,x,"3-path") for x in findall(==(1),rel)])
    #jpaths
    append!(ships,[(0,i,j,x,"3-path") for x in findall(==(2),rel)])
    #triangles
    append!(ships,[(0,i,j,x,"3-tri") for x in findall(==(3),rel)])

    #4-paths iedge
    append!(ships,[(j,i,x,y,"4-path") for x in findall(==(1),rel) for y in last.(filter(z->first(z)==x,Rel[4]))])
    #4-paths jedge
    append!(ships,[(i,j,x,y,"4-path") for x in findall(==(2),rel) for y in last.(filter(z->first(z)==x,Rel[4]))])



    #4-tails icentre
    append!(ships,[(y,x,i,j,"4-tail") for x in findall(==(1),rel) for y in last.(filter(z->first(z)==x,Rel[6]))])
    #4-tails jcentre
    append!(ships,[(y,x,j,i,"4-tail") for x in findall(==(2),rel) for y in last.(filter(z->first(z)==x,Rel[6]))])
    #4-tails tricentre
    append!(ships,[(i,j,x,y,"4-tail") for x in findall(==(3),rel) for y in last.(filter(z->first(z)==x,Rel[6]))])

    #4-cycles
    append!(ships,[(i,j,x,y,"4-cycle") for x in findall(==(2),rel) for y in last.(filter(z->first(z)==x,Rel[7]))])
    
    #4-chord iedge orbit
    append!(ships,[(j,i,x,y,"4-chord") for x in findall(==(3),rel) for y in last.(filter(z->first(z)==x,Rel[8]))])
    #4-chord jedge orbit
    append!(ships,[(i,j,x,y,"4-chord") for x in findall(==(3),rel) for y in last.(filter(z->first(z)==x,Rel[9]))])

    #4-clique
    append!(ships,[(x,i,j,y,"4-clique") for x in findall(==(3),rel) for y in last.(filter(z->first(z)==x,Rel[10]))])

    ##combinatorials might be more difficult...
    #4-paths centre orbit
    append!(ships,[(x,i,j,y,"4-path") for x in findall(==(1),rel) for y in findall(==(2),rel) if !((y,x) in Rel[7])])
    
    #4-chords centre orbit
    append!(ships,[(x,i,j,y,"4-chord") for x in findall(==(3),rel) for y in findall(==(3),rel) if !(y<x && (x=>y) in Rel[10])])
    
    #4-stars i centre
    append!(ships,[(x,j,i,y,"4-star") for x in findall(==(1),rel) for y in findall(==(1),rel) if !(y<x && (x=>y) in Rel[6])])
    
    #4-stars j centre
    append!(ships,[(x,i,j,y,"4-star") for x in findall(==(2),rel) for y in findall(==(2),rel) if !(y<x && (x=>y) in Rel[6])])

    #4-tails tri i edge
    append!(ships,[(x,j,i,y,"4-tail") for x in findall(==(3),rel) for y in findall(==(1),rel) if !((x=>y) in Rel[8])])
    
    #4-tails tri j edge
    append!(ships,[(x,i,j,y,"4-tail") for x in findall(==(3),rel) for y in findall(==(2),rel) if !((x=>y) in Rel[9])])
end
    return ships

end

function per_edge_counts(edge::Int,vertex_type_list::Vector{<:AbstractString},edgelist::Union{Array{Pair{Int,Int},1},Array{Pair,1}},graphlet_size::Int,neighbourdict::Dict{Int,Vector{Int}})
    count_dict = DefaultDict{String,Int}(0)
    h=edge  
#   # get nodes i and j for this edge
        i = edgelist[h].first
        j = edgelist[h].second
        #get neighbourhoods of i and j
        gamma_i = neighbourdict[i]  
        gamma_j = neighbourdict[j]
        
        delim = "_"
        #three node graphlets
    ##more efficient loop based method
    #rel = DefaultDict{Int,Int}(0)
    rel = zeros(Int,length(vertex_type_list))
    Tri = Array{Int,1}()
    iPath = Array{Int,1}()
    jPath = Array{Int,1}()
    #prelim cycle through i neighbours
    for k in gamma_i
        if (k!=j)
            rel[k] = 1
        end
    end
    for k in gamma_j
        if(k!=i)
            if (rel[k]==1)
                ##triangle
                count_dict[graphlet_string(vertex_type_list[i],vertex_type_list[j],vertex_type_list[k],"3-tri",delim)]+=1
                append!(Tri,k)
                rel[k] = 3
            else
                #j-path
                count_dict[graphlet_string(vertex_type_list[i],vertex_type_list[j],vertex_type_list[k],"3-path",delim)]+=1
                append!(jPath,k)
                rel[k] = 2
            end
        end
    end
        
    for k in gamma_i
        if (rel[k]==1)
            #ipaths
            count_dict[graphlet_string(vertex_type_list[j],vertex_type_list[i],vertex_type_list[k],"3-path",delim)]+=1
            append!(iPath,k)
        end
    end 
        if (graphlet_size==4)
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
                        count_dict[graphlet_string(a,b,vertex_type_list[i],vertex_type_list[j],"4-star",delim)] += Int(0.5*iPathTypes[inda]*(iPathTypes[inda]-1)) -     count_dict[graphlet_string(a,b,vertex_type_list[i],vertex_type_list[j],"4-tail-edge-orbit",delim)]
                        #j-centre
                        count_dict[graphlet_string(a,b,vertex_type_list[j],vertex_type_list[i],"4-star",delim)] += Int(0.5*jPathTypes[inda]*(jPathTypes[inda]-1)) -     count_dict[graphlet_string(a,b,vertex_type_list[j],vertex_type_list[i],"4-tail-edge-orbit",delim)]
                        ##4-tail tri-edge orbit
                        #i-centre
                        count_dict[graphlet_string(a,vertex_type_list[j],vertex_type_list[i],b,"4-tail-tri-edge-orbit",delim)] += TriTypes[inda]*iPathTypes[inda] -     count_dict[graphlet_string(a,vertex_type_list[i],b,vertex_type_list[j],"4-chord-edge-orbit",delim)]
                        #j-centre
                        count_dict[graphlet_string(a,vertex_type_list[i],vertex_type_list[j],b,"4-tail-tri-edge-orbit",delim)] += TriTypes[inda]*jPathTypes[inda] -     count_dict[graphlet_string(a,vertex_type_list[j],b,vertex_type_list[i],"4-chord-edge-orbit",delim)]

                    else # when i and j are also of same type
                        ##4-star
                        #To maintain type order here, we also have to separate. Note we enforce that the centre of the star is THIRD listed (in line with the 4-tail edge orbit layout)
                        count_dict[graphlet_string(a,b,vertex_type_list[i],vertex_type_list[j],"4-star",delim)] += Int(0.5*iPathTypes[inda]*(iPathTypes[inda]-1))+Int(0.5*jPathTypes[inda]*(jPathTypes[inda]-1)) -  count_dict[graphlet_string(a,b,vertex_type_list[i],vertex_type_list[j],"4-tail-edge-orbit",delim)]
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
                    else ##i and j are same type
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

#aggregator function
function t2(d1,d2)
    append!(d1,d2)
    d1
end

function count_graphlets(args...;run_method::String="serial",progress::Bool=false)
    Chi = local_graphlets(args...;run_method=run_method,progress=progress)
    

    graphlets = total_graphlets(Chi)
    return graphlets
end


function graphlet_relationships(vertex_type_list::Vector{<:AbstractString},edgelist::Union{Array{Pair{Int,Int},1},Array{Pair,1}},graphlet_size::Int=3;run_method::String="serial",temp_dir::String="rel_dir",progress::Bool=false)
    #get neighbourhood for each vertex in advance (rather than calling per-edge)
    neighbourdict=Neighbours(edgelist)
    #preallocate array to store each edge relationship dict 
    #Rel=Array{Array{Tuple{Int,Int,Int,Int,String},1}}(undef,size(edgelist,1));
    #Rel = Array{Int,3}(undef,length(vertex_type_list),length(vertex_type_list),length(vertex_type_list))

    ##for larger/more connected networks, counting relationships will only be possible if outputs are written directly to a CSV file. This is now enforced for all relationship counting
    #remove any existing temp dir and remake
    run(`rm -rf $temp_dir`)
    run(`mkdir $temp_dir`)
    #now run per edge relationship counts, generating relationship csv files in temp_dir (one for each worker process)
    if(run_method == "distributed")
        if(progress == true)
            @info "Distributing edges to workers"

            @showprogress pmap(x->per_edge_relationships(x,vertex_type_list,edgelist,graphlet_size,neighbourdict,temp_dir,write_out=true),1:size(edgelist,1),batch_size =1000)
        else
            pmap(x->per_edge_relationships(x,vertex_type_list,edgelist,graphlet_size,neighbourdict,temp_dir,write_out=true),1:size(edgelist,1),batch_size =1000)        
        end
    elseif (run_method == "serial")
        for h in 1 :size(edgelist,1)
            per_edge_relationships(h,vertex_type_list,edgelist,graphlet_size,neighbourdict,temp_dir,write_out=true)
        end
    else
        throw(ArgumentError("run_method not recognised."))
    end
    
    ## we now need to collate and sort these relationship files, removing duplicates. This still needs to be done on disk memory

    ##first sort into graphlet type so that all duplicates are guaranteed to be in same file
    progress ? (@info "Sorting into graphlet types for deduplication") : nothing 
    for file  in filter(x->occursin(".csv",x),readdir(temp_dir,join=true))
        #sort file into each graphlet (appending to existing if it exists)
        run(`awk -v var="$(temp_dir)/" -F, '{print >> (var $NF".csv")}' $file`)
        #remove original file (for space constraints)
        run(`rm $file`)
    end

    ##now sort each graphlet file, based on the orbit structure of that graphlet
   #TODO update to use perl? and perhaps improve parallelisation  
    progress ? (@info "Sorting each graphlet type according to orbit symmetries") : nothing 
     @sync @distributed for file in filter(x->occursin(".csv",x),readdir(temp_dir,join=true))
        if(splitext(basename(file))[1] == "3-path")
            ## remove 0 field and sort first and third nodes
            run(`awk -v var="$(temp_dir)/" -F, '$2<$4{print $1,$2,$3,$4,$5 > var $5"_sorted.csv"} $4<$2{print $1,$4,$3,$2,$5 > var $5"_sorted.csv"}' $file`) 

        elseif(splitext(basename(file))[1] == "3-tri")
            ## remove 0 field and sort all three nodes
            run(`awk -v var="$(temp_dir)/" -F, '
                $2<$3 && $2<$4 && $3<$4{print $1,$2,$3,$4,$5 > var $5"_sorted.csv"}
                $2<$3 && $2<$4 && $4<$3{print $1,$2,$4,$3,$5 > var $5"_sorted.csv"}
                $3<$2 && $2<$4 && $3<$4{print $1,$3,$2,$4,$5 > var $5"_sorted.csv"}
                $3<$2 && $4<$2 && $3<$4{print $1,$3,$4,$2,$5 > var $5"_sorted.csv"}
                $2<$3 && $4<$2 && $4<$3{print $1,$4,$2,$3,$5 > var $5"_sorted.csv"}
                $3<$2 && $4<$2 && $4<$3{print $1,$4,$3,$2,$5 > var $5"_sorted.csv"}
                ' $file`) 
        elseif(splitext(basename(file))[1] == "4-path")
            ## sort inner, and then sort outer accordingly
            run(`awk -v var="$(temp_dir)/" -F, '
                $2<$3{print $1,$2,$3,$4,$5 > var $5"_sorted.csv"}
                $3<$2{print $4,$3,$2,$1,$5 > var $5"_sorted.csv"}
                ' $file`) 
        elseif(splitext(basename(file))[1] == "4-star")
            ## sort on 1,2 and 4
            run(`awk -v var="$(temp_dir)/" -F, '
                $1<$2 && $1<$4 && $2<$4{print$1,$2,$3,$4,$5 > var $5"_sorted.csv"}
                $1<$2 && $1<$4 && $4<$2{print$1,$4,$3,$2,$5 > var $5"_sorted.csv"}
                $2<$1 && $1<$4 && $2<$4{print$2,$1,$3,$4,$5 > var $5"_sorted.csv"}
                $2<$1 && $4<$2 && $2<$4{print$2,$4,$3,$1,$5 > var $5"_sorted.csv"}
                $1<$2 && $4<$1 && $4<$2{print$4,$1,$3,$2,$5 > var $5"_sorted.csv"}
                $1<$2 && $1<$4 && $2<$4{print$4,$2,$3,$1,$5 > var $5"_sorted.csv"}
                ' $file`) 
        elseif(splitext(basename(file))[1] == "4-tail")
            ## sort on 1 and 2 only
                run(`awk -v var="$(temp_dir)/" -F, '
                $1<$2{print $1,$2,$3,$4,$5 > var $5"_sorted.csv"}
                $2<$1{print $2,$1,$3,$4,$5 > var $5"_sorted.csv"}
                ' $file`) 
        elseif(splitext(basename(file))[1] == "4-cycle")
            ## sort on all to orientate on lowest (in pos 1), then sort adjacents for position 2 and 4, and then position node non-adjacent to position 3
                run(`awk -v var="$(temp_dir)/" -F, '
                $1<$2 && $1<$3 && $1<$4 && $2<$4  {print$1,$2,$3,$4,$5 > var $5"_sorted.csv"}
                $1<$2 && $1<$3 && $1<$4 && $4<$2  {print$1,$4,$3,$2,$5 > var $5"_sorted.csv"}
                $2<$1 && $2<$3 && $2<$4 && $1<$3  {print$2,$1,$4,$3,$5 > var $5"_sorted.csv"}
                $2<$1 && $2<$3 && $2<$4 && $3<$1  {print$2,$3,$4,$1,$5 > var $5"_sorted.csv"}
                $3<$1 && $3<$2 && $3<$4 && $2<$4  {print$3,$2,$1,$4,$5 > var $5"_sorted.csv"}
                $3<$1 && $3<$2 && $3<$4 && $4<$2  {print$3,$4,$1,$2,$5 > var $5"_sorted.csv"}
                $4<$1 && $4<$2 && $4<$3 && $1<$3  {print$4,$1,$2,$3,$5 > var $5"_sorted.csv"}
                $4<$1 && $4<$2 && $4<$3 && $3<$1  {print$4,$3,$2,$1,$5 > var $5"_sorted.csv"}
                ' $file`) 
        elseif(splitext(basename(file))[1] == "4-chord")
            ## sort inner and outer independently
            run(`awk -v var="$(temp_dir)/" -F, '
                $2<$3 && $1<$4{print $1,$2,$3,$4,$5 > var $5"_sorted.csv"}
                $2<$3 && $4<$1{print $4,$2,$3,$1,$5 > var $5"_sorted.csv"}
                $3<$2 && $1<$4{print $1,$3,$2,$4,$5 > var $5"_sorted.csv"}
                $3<$2 && $4<$1{print $4,$3,$2,$1,$5 > var $5"_sorted.csv"}
                ' $file`) 
        elseif(splitext(basename(file))[1] == "4-clique")
            ## sort all nodes
            run(`awk -v var="$(temp_dir)/" -F, '
                $1<$2 && $1<$3 && $1<$4 && $2<$3 && $2<$4 && $3<$4  {print $1,$2,$3,$4,$5 > var $5"_sorted.csv"}
                $1<$2 && $1<$3 && $1<$4 && $2<$3 && $2<$4 && $4<$3  {print $1,$2,$4,$3,$5 > var $5"_sorted.csv"}
                $1<$2 && $1<$3 && $1<$4 && $3<$2 && $3<$4 && $2<$4  {print $1,$3,$2,$4,$5 > var $5"_sorted.csv"}
                $1<$2 && $1<$3 && $1<$4 && $3<$2 && $3<$4 && $4<$2  {print $1,$3,$4,$2,$5 > var $5"_sorted.csv"}
                $1<$2 && $1<$3 && $1<$4 && $4<$2 && $4<$3 && $2<$3  {print $1,$4,$2,$3,$5 > var $5"_sorted.csv"}
                $1<$2 && $1<$3 && $1<$4 && $4<$2 && $4<$3 && $3<$2  {print $1,$4,$3,$2,$5 > var $5"_sorted.csv"}
                $2<$1 && $2<$3 && $2<$4 && $1<$3 && $1<$4 && $3<$4  {print $2,$1,$3,$4,$5 > var $5"_sorted.csv"}
                $2<$1 && $2<$3 && $2<$4 && $1<$3 && $1<$4 && $4<$3  {print $2,$1,$4,$3,$5 > var $5"_sorted.csv"}
                $2<$1 && $2<$3 && $2<$4 && $3<$1 && $3<$4 && $1<$4  {print $2,$3,$1,$4,$5 > var $5"_sorted.csv"}
                $2<$1 && $2<$3 && $2<$4 && $3<$1 && $3<$4 && $4<$1  {print $2,$3,$4,$1,$5 > var $5"_sorted.csv"}
                $2<$1 && $2<$3 && $2<$4 && $4<$1 && $4<$3 && $1<$3  {print $2,$4,$1,$3,$5 > var $5"_sorted.csv"}
                $2<$1 && $2<$3 && $2<$4 && $4<$1 && $4<$3 && $3<$1  {print $2,$4,$3,$1,$5 > var $5"_sorted.csv"}
                $3<$1 && $3<$2 && $3<$4 && $1<$2 && $1<$4 && $2<$4  {print $3,$1,$2,$4,$5 > var $5"_sorted.csv"}
                $3<$1 && $3<$2 && $3<$4 && $1<$2 && $1<$4 && $4<$2  {print $3,$1,$4,$2,$5 > var $5"_sorted.csv"}
                $3<$1 && $3<$2 && $3<$4 && $2<$1 && $2<$4 && $1<$4  {print $3,$2,$1,$4,$5 > var $5"_sorted.csv"}
                $3<$1 && $3<$2 && $3<$4 && $2<$1 && $2<$4 && $4<$1  {print $3,$2,$4,$1,$5 > var $5"_sorted.csv"}
                $3<$1 && $3<$2 && $3<$4 && $4<$1 && $4<$2 && $1<$2  {print $3,$4,$1,$2,$5 > var $5"_sorted.csv"}
                $3<$1 && $3<$2 && $3<$4 && $4<$1 && $4<$2 && $2<$1  {print $3,$4,$2,$1,$5 > var $5"_sorted.csv"}
                $4<$1 && $4<$2 && $4<$3 && $1<$2 && $1<$3 && $2<$3  {print $4,$1,$2,$3,$5 > var $5"_sorted.csv"}
                $4<$1 && $4<$2 && $4<$3 && $1<$2 && $1<$3 && $3<$2  {print $4,$1,$3,$2,$5 > var $5"_sorted.csv"}
                $4<$1 && $4<$2 && $4<$3 && $2<$1 && $2<$3 && $1<$3  {print $4,$2,$1,$3,$5 > var $5"_sorted.csv"}
                $4<$1 && $4<$2 && $4<$3 && $2<$1 && $2<$3 && $3<$1  {print $4,$2,$3,$1,$5 > var $5"_sorted.csv"}
                $4<$1 && $4<$2 && $4<$3 && $3<$1 && $3<$2 && $1<$2  {print $4,$3,$1,$2,$5 > var $5"_sorted.csv"}
                $4<$1 && $4<$2 && $4<$3 && $3<$1 && $3<$2 && $2<$1  {print $4,$3,$2,$1,$5 > var $5"_sorted.csv"}
                ' $file`) 
        else
            throw(ArgumentError("Symmetries have not yet been defined for this graphlet"))
        end    
        #remove original file (for space constraints)
        run(`rm $file`)
    end

    ##now remove duplicates from each graphlet file, writing to one relationships file
    #TODO this is currently breaking on larger networks, where one graphlet with lots of occurences will break memory constraints
    #solution is to break up the deduplication into batches, and then condensing into one smaller file, where final (cross bin) duplicates can be found (is this best method though?)
    progress ? (@info "Removing duplicates") : nothing 
    for file in filter(x->occursin(".csv",x),readdir(temp_dir,join=true))
        #run(`/bin/bash -c "sort -u $(file) >> $(temp_dir)/relationships.csv"`)
        #run(`awk -v var="$(temp_dir)/" '!seen[$0] {print >> var"relationships.csv"} {seen[$0] += 1}' $file`)
        ### for sort: if switch to individual output files use -o $(temp_dir)/$(split(file,".")[1])_dedup.csv
        #run(`sort -us --parallel=$(Threads.nthreads()) -T $(temp_dir) $file '>>' $(temp_dir)/relationships.csv`)
        open( io -> run(pipeline(`sort -us --parallel=$(Threads.nthreads()) -T $(temp_dir) $file`,stdout=io)), "$(temp_dir)/relationships.csv","a")
        run(`rm $file`)
    end
    
    @info "Finished. Relationships stored in a CSV file at $(temp_dir)/relationships.csv"    

    ## load final file onto heap memory
    
    #progress ? (@info "Loading relationships into memory") : nothing 
    #Rel = CSV.read("$temp_dir/relationships.csv",DataFrame,header=false)
    return nothing 
end


function local_graphlets(vertex_type_list::Vector{<:AbstractString},edgelist::Union{Array{Pair{Int,Int},1},Array{Pair,1}},graphlet_size::Int=3;run_method::String="serial",progress::Bool=false)
    ##INPUTS TO PER EDGE FUNCTION
    #get neighbourhood for each vertex in advance (rather than calling per-edge)
    neighbourdict=Neighbours(edgelist)
    # set up function to apply dict to arrays (maybe a better way exists, but this works?) DO NOT NEED AT PRESENT
    #neighbourdictfunc(x::Int) = neighbourdict[x]
    #vertex set derived from edgelist (allows for any arbritary vertex labelling)
    
    #preallocate array to store each edge's graphlet dictionary 
    Chi=Array{Dict{String,Int}}(undef,size(edgelist,1));
    #per edge process
    if(run_method == "threads")
        throw(ArgumentError("threads method is currently not working."))
        ##NOT WORKING AT PRESENT
        Threads.@threads for h in 1 :size(edgelist,1)
            edge = per_edge_counts(h,vertex_type_list,edgelist,graphlet_size,neighbourdict)
            Chi[h] = edge[1]
        end
    elseif(run_method == "remote-channel")
        
        ##run tasks through input and ooutput remote-channels (r is output remote channel)
        r = remote_channel_method(vertex_type_list,edgelist,graphlet_size,neighbourdict;progress = progress)
        ## there should be length(edgelist) results on channel r. we harvest them into each array as required:
        if(progress)
            @info "Taking results from output channel..."
            @showprogress for n in 1:length(edgelist)
                i,chi = take!(r)
                Chi[i] = chi
            end
        else
            for n in 1:length(edgelist)
                i,chi = take!(r)
                Chi[i] = chi
            end
        end

    elseif(run_method == "distributed")
        if(progress == true)

            @info "Distributing edges to workers..."
            ##alternative option using pmap (dynamically manages worker loads, so that all CPUS are used for entire job. Needs some mechanism for reduction at end though
            Chi = @showprogress pmap(x->per_edge_counts(x,vertex_type_list,edgelist,graphlet_size,neighbourdict),1:size(edgelist,1),batch_size =1000)
        else
            Chi = pmap(x->per_edge_counts(x,vertex_type_list,edgelist,graphlet_size,neighbourdict),1:size(edgelist,1),batch_size =1000)
        end
    elseif(run_method == "distributed-old")
        if(progress==true)
            @info "Distributing edges to workers..."

            res = @showprogress @distributed (t2) for h in 1:size(edgelist,1)
                [(h,per_edge_counts(h,vertex_type_list,edgelist,graphlet_size,neighbourdict))]        
            end
        else
            res = @sync @distributed (t2) for h in 1:size(edgelist,1)
                [(h,per_edge_counts(h,vertex_type_list,edgelist,graphlet_size,neighbourdict))]        
            end
        end
        for r in res
            Chi[first(r)] = last(r)
        end

    elseif (run_method == "serial")
        for h in 1 :size(edgelist,1)
            edge = per_edge_counts(h,vertex_type_list,edgelist,graphlet_size,neighbourdict)
            Chi[h] = edge
        end
    else
        throw(ArgumentError("run_method not recognised."))
    end
    return Chi
end

function total_graphlets(Chi::Array{AbstractDict{String,Int}})    
    @info "Calculating total counts for each graphlet..."
    #total counts for each graphlet
    total_counts = reduce(mergecum,Chi)

    #reorder names to merge orbits
    graphlet_names = (split.(collect(keys(total_counts)),"_"))
    for el in 1:size(graphlet_names,1)
        #for 3 graphlets:
        if(length(graphlet_names[el])==4)

            #for x-y-z paths such that y!=x AND y!=z (different orbit to other 3-paths) we reorder only the end nodes
            if (graphlet_names[el][1]!=graphlet_names[el][2] && graphlet_names[el][3]!=graphlet_names[el][2] && graphlet_names[el][4]=="3-path")
                graphlet_names[el][[1,3]] =sort(graphlet_names[el][[1,3]])
            else
                #for all other paths and triangles
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
                #we iwant to switch if exterior needs switching: 
                if (graphlet_names[el][[1,4]][1] != sort(graphlet_names[el][[1,4]])[1])
                    graphlet_names[el][[1,4]] = sort(graphlet_names[el][[1,4]])
                    graphlet_names[el][[2,3]] = graphlet_names[el][[3,2]]
                end
                #or if exterior is the same, we sort inner types:
                if (graphlet_names[el][[1]] == graphlet_names[el][[4]])
                    graphlet_names[el][[2,3]] = sort(graphlet_names[el][[2,3]]) 
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
                ## cycles (differentiate between instances where a) there is a triplet of same type, which will occupy first three slots; b) at least on pair of adjacent nodes are of same type, in which case make first two be the pair (and if there are two pairs, choose first pair via sorting),; or c) all adjacent node pairs are of different type, in which case let any matching non adjacent pairs take slots 1 and 3, and if there two non adjacent pairs, choose first pair via sorting).    
            elseif (graphlet_names[el][5] == "4-cycle")
                #get occurences of each type in graphlet
                occs = countmap(graphlet_names[el][1:4])
                moccs = max(values(occs)...)
                if (moccs == 4)
                    # no switching needed (homogeneous graphlet) 
                elseif (moccs == 3)
                    #triplet case
                    #find type and position that is not in triplet
                    solo = collect(keys(filter(x->last(x) !== moccs,occs)))[1]
                    switch = findall(x->x == solo,graphlet_names[el])
                    #switch solo type to last position (potentially trivially but thats ok) 
                    graphlet_names[el][switch...] = graphlet_names[el][4]
                    graphlet_names[el][4] = solo
                elseif (moccs == 2)
                    ## find which first type that matches max occs 
                    primary_pair_type = sort(collect(keys(filter(x->last(x) == moccs,occs))))[1]
                    sig = graphlet_names[el].==primary_pair_type
                    ##check for case c)
                    if (sig == BitArray(vec([true false true false false])))
                        ## sort non primary non-adjecent pair nodes 
                        graphlet_names[el][[2,4]] = sort(graphlet_names[el][[2,4]])
                    elseif (sig == BitArray(vec([false true false true false])))
                        ##switch primary nonadjacent pair to positions 1 and 3, and sort the other nodes
                        graphlet_names[el][[2,4]] = sort(graphlet_names[el][[1,3]])
                        graphlet_names[el][[1,3]] = [primary_pair_type, primary_pair_type]
                    else
                        ## it must be case b)
                        graphlet_names[el][[3,4]] = sort(filter(x->x!=primary_pair_type,graphlet_names[el][1:4]))
                        graphlet_names[el][[1,2]] = [primary_pair_type, primary_pair_type]

                    end
                else
                    ## all types must be unique (and thus there must be at least 4 types present) so we sort on all nodes
                    graphlet_names[el][1:4] = sort(graphlet_names[el][1:4])

                    #now need to find out if 
                    #if (length(occs) == 2)

                end
                #for cliques, just order everything
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

function find_motifs(edgelist::Union{Array{Pair{Int,Int},1},Array{Pair,1}},null_model::String,null_num::Int; typed::Bool=false, typelist::Array{String,1}=nothing,plotfile::String="DONOTPLOT",graphlet_size::Int=3)
    ##Calculate null model counts
    edgelists = edgelists_null_model(edgelist,null_num,null_model,typelist,graphlet_size)
    null_model_calc = null_model_counts(typelist,edgelists)
    null_model_df = null_model_dataframe(null_model_calc) 
    
    # calculate real network counts
    graphlet_counts = count_graphlets(typelist,edgelist,graphlet_size,run_method="distributed")[1]
    
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

        #get graphlet counts in df form
        graphlet_df = DataFrame(graphlet = broadcast(first,collect(graphlet_counts)),value = broadcast(last,collect(graphlet_counts)))
        #convert to concentrations for plotting
        null_model_conc = null_model_concentrations(null_model_calc)
        conc_null_model_df = null_model_dataframe(null_model_conc) 
        graphlet_concentrations = concentrate(graphlet_counts) 
        conc_graphlet_df = DataFrame(graphlet = broadcast(first,collect(graphlet_concentrations)),value = broadcast(last,collect(graphlet_concentrations)))
        
        
        #trying log version (using counts not concentrations)
        log_null_model_df = copy(null_model_df)
        log_null_model_df.value = log.(null_model_df.value)
        log_graphlet_df = copy(graphlet_df)
        log_graphlet_df.value = log.(graphlet_df.value)
        #trying log log version (using counts not concentrations)
        log_log_null_model_df = copy(log_null_model_df)
        log_log_null_model_df.value = log.(log_null_model_df.value)
        log_log_graphlet_df = copy(log_graphlet_df)
        log_log_graphlet_df.value = log.(graphlet_df.value)
        #setting up plot matrix
        hom_graphlets = unique(last.(split.(graphlet_df[:graphlet],"_")))
        xdim = 3
        ydim = Int(ceil(length(hom_graphlets)/xdim))
        plots=Array{Union{Plot,Context}}(undef,xdim,ydim)
        for (i,g) in enumerate(hom_graphlets)
            ##without graphlet points (test)
        #   plots[i] = plot(layer(filter(:graphlet=>x->occursin("coding_coding_coding_"*g,x),filter(:graphlet=>x->occursin(g,x),null_model_df)),x=:graphlet,y=:value,Geom.boxplot(suppress_outliers = false),color=:graphlet),Guide.xticks(label=true),Theme(key_position = :none));
        plots[i] = plot(layer(filter(:graphlet=>x->occursin(g,x),log_null_model_df),x=:graphlet,y=:value,Geom.boxplot(suppress_outliers = true),color=:graphlet),layer(filter(:graphlet=>x->occursin(g,x),log_graphlet_df),x = :graphlet,y = :value, Geom.point,color=["count in graph"]),Guide.xticks(label=true),Theme(key_position = :none),Guide.xlabel(nothing),Guide.ylabel("log value"),Guide.yticks(orientation=:vertical));
        end
        #add empty panels
        for i in length(hom_graphlets)+1:length(plots)
            plots[i] = context()
        end
        draw(SVG(plotfile,10inch,20inch),gridstack(plots))
        @info "Stat plots saved to $plotfile."
    end
    #null_model_sum = reduce(mergecum,null_model_calc)
    return [zscores,edgelists,null_model_calc,graphlet_counts,plots]
end 

#put function everywhere
function no_work(jobs, results,args...) # define work function everywhere
    while true
        job_ids = take!(jobs)
        for i in job_ids
            Chi = GraphletCounting.per_edge_counts(i,args...)
            put!(results,(i,Chi))
        end
     end
end

function do_work(jobs, results,args...) # define work function everywhere
    while true
        job_id = take!(jobs)
        Chi = GraphletCounting.per_edge_counts(job_id,args...)
        put!(results,(job_id,Chi))
    end
end

function remote_channel_method(args...;progress::Bool=false)
    
    batch_size = 1
    m = length(args[2])
    #create input and output channels
    #jobs = RemoteChannel(()->Channel{Vector{Int}}(m));
    jobs = RemoteChannel(()->Channel{Int}(m));
    results = RemoteChannel(()->Channel{Tuple}(m));
    #add jobs to input channel
    for i in 1:batch_size:m
        #put!(jobs,collect(i:(i+batch_size-1)))
        put!(jobs,i)
    end
    if (progress)
        @info "Distributing edges to workers..."
        @showprogress for p in workers() # start tasks on the workers to process requests in parallel
            remote_do(GraphletCounting.no_work, p, jobs, results,args...)
        end
    else
        for p in workers() # start tasks on the workers to process requests in parallel
            remote_do(GraphletCounting.do_work, p, jobs, results,args...)
        end

    end
return results
end
