using LightGraphs,RCall,DataFrames,GraphletCounting,NetworkConstruction
#using PrettyTables

function connected_components_html_table(adjacency_matrix::AbstractArray,filename::String)
    g = Graph(adjacency_matrix)
    cc = size.(connected_components(g),1)
    io = open(filename, "w")
    println(io,"""<!DOCTYPE html>
<html>
<meta charset="UTF-8">
<style>
table, td, th {
    border-collapse: collapse;
    font-family: sans-serif;
}

td, th {
    border-bottom: 0;
    padding: 4px
}

tr:nth-child(odd) {
    background: #eee;
}

tr:nth-child(even) {
    background: #fff;
}

tr.header {
    background: navy !important;
    color: white;
    font-weight: bold;
}

tr.subheader {
    background: lightgray !important;
    color: black;
}

tr.headerLastRow {
    border-bottom: 2px solid black;
}

th.rowNumber, td.rowNumber {
    text-align: right;
}

</style>
<body>""")
    println(io,"""<p style="margin-bottom:3cm;">
        <big> Number of vertices: """,size(adjacency_matrix,1),"""</big>
        </p>""")
    println(io,"""<p style="margin-bottom:3cm;">
        <big> Number of edges: """,Int(sum(adjacency_matrix)/2),"""</big>
        </p>""")
    println(io,"""<p style="margin-bottom:3cm;">
        <big> Connected Components:</big>
        </p>""")
    header = string.(1:length(cc))
    #table = pretty_table(String,cc', vec(header), backend = :html,standalone = false);
    println(io,table)
    println(io,"""</body>
        </html>""")
    close(io)
end

function get_community_structure(adj_matrix::AbstractArray,vertex_names::Array{String,1},detection_method::String;plot_prefix::String = "", threejs_plot::Bool = false)
    ## Get community structure of a network using r-igraph package 
    R"""
    sapply(names(sessionInfo()$otherPkgs),function(pkg) detach(paste0('package:',pkg),character.only =T,force = T));rm(list=ls())
    """
    @rput adj_matrix
    @rput vertex_names
    @rput detection_method
    @rput threejs_plot
    @rput plot_prefix

    R"""
    library(RColorBrewer)
    library(igraph)
    library(threejs)
    library(htmlwidgets)
    #keep this last to mask other packages:
    library(tidyverse)
    
    ##change to tibble
    vertex_names = tibble(vertex_names)
    vertex_names_trimmed = sapply(vertex_names,tools::file_path_sans_ext)
    pre_g <- graph.adjacency(adj_matrix,mode="undirected",diag=F)
    #Now we add attributes to graph by rebuilding it via edge and vertex lists
    edges <- as_tibble(as.data.frame(get.edgelist(pre_g)))
    vertices <- tibble(name = 1:nrow(vertex_names),vertex_names)%>% dplyr::rename(label=vertex_names)
    
    g <- graph_from_data_frame(edges,directed = F, vertices)
    #
    ##delete zero degree vertices
    g <- delete.vertices(g,igraph::degree(g)==0)
    ##find maximum connected component
    g <- decompose(g, mode = "strong", max.comps = NA, min.vertices = 10)[[1]]
    #
    ###Community structure
    ##get only connected vertices
    vertices <- vertices %>% filter(name %in% names(V(g)))
    
    ##cluster graph and colour by community
    if(detection_method == "louvain")
            {
                communities <- cluster_louvain(g)
            }
    
    ##add colours to graph
    vertices <- vertices %>% mutate(group = as_factor(communities$membership),color = group)
    ##if there are more than 11 communities, spectral colour palette is not sufficient. so we concat two palettes (if there are more than 22 comms, need another!)
    n =  length(unique(communities$membership)) 
    if(n<12)
    {
        colour_palette = brewer.pal(name = "Spectral", n = n)
    } else if (n>12 & n<23)
    {
        x =n-11
        colour_palette = c(brewer.pal(name = "Spectral", n = 11),brewer.pal(name = "BrBG", n = x))
    }else if (n>23 & n<34)
    {
        x =n-22
        colour_palette = c(brewer.pal(name = "Spectral", n = 11),brewer.pal(name = "BrBG", n = 11),brewer.pal(name = "RdYlBu", n = x))
    }else if (n>34 & n<45)
    {
        x =n-33
        colour_palette = c(brewer.pal(name = "Spectral", n = 11),brewer.pal(name = "BrBG", n = 11),brewer.pal(name = "RdYlBu", n = 11),brewer.pal(name="PiYG",n=x))
    }else
    {
        paste("ERROR: too many communiies to colour network")
    }

    levels(vertices$color) <- colour_palette
    edges <- as_tibble(as.data.frame(get.edgelist(g)))
    g <- graph_from_data_frame(edges,directed = F, vertices)
    vertex_attr(g,"size") <- 0.5
    
    if(threejs_plot == T)
        {
    
            plot = graphjs(g,bg = "white");
            saveWidget(plot,paste0(plot_prefix,"_communities.html"))
        }
    """
    
    @rget vertices
    
    return vertices
end

function get_community_node_types(adj_matrix::AbstractArray,community_vector::Array{Int,1},type_vector::Array{String,1})
    if (size(adj_matrix,1) != length(community_vector))
        @error "Community vector is of different size to adjacency matrix. Some nodes do not have a community assigned (use 0 for nodes in no community)."
        end

    #modify adjacency matrix so that every non-zero entry indicates the community of the connecting (columnwise) node and every zero entry indicates the community of the (rowwise) node
    comm_matrix = (adj_matrix.*community_vector)'+(.!adj_matrix.*community_vector)
    ##the communities that each nodes is connected to
    connections = unique.(eachrow(comm_matrix))
    ## all nodes that are only connected to their own community
    interior_nodes = length.(connections).==1
    #per community:
    comm_df = DataFrame()
    for comm in unique(community_vector)
        ##boolean of nodes in community
        comm_nodes = community_vector.== comm
        ## nodes that are connected to nodes in that community.
        conn_nodes = in.(comm,connections)
        ##nodes in commmunity that are only connected community
        comm_interiors = comm_nodes.*interior_nodes
        ##nodes in community that are connected to nodes in another community
        comm_boundaries = comm_nodes.*.!interior_nodes
        ##nodes not in community that are connected to nodes in community
        comm_neighbours = conn_nodes.*.!comm_nodes
        ##nodes not in community that are not connected to nodes in community
        comm_exteriors = .!conn_nodes.*.!comm_nodes
        append!(comm_df,DataFrame(Community = comm, Size = sum(comm_nodes),Coding_Interiors = sum(type_vector[comm_interiors].=="coding" ),Noncoding_Interiors = sum(type_vector[comm_interiors].=="noncoding" ),Coding_Boundaries = sum(type_vector[comm_boundaries].=="coding" ),Noncoding_Boundaries = sum(type_vector[comm_boundaries].=="noncoding" ),Coding_Neighbours = sum(type_vector[comm_neighbours].=="coding" ),Noncoding_Neighbours = sum(type_vector[comm_neighbours].=="noncoding" ),Coding_Exteriors = sum(type_vector[comm_exteriors].=="coding" ),Noncoding_Exteriors = sum(type_vector[comm_exteriors].=="noncoding" )))

    end
    return comm_df
end

function get_functional_annotations(comm_vertices::DataFrame;ensembl_version::String="current",write_csv::Bool = false,csv_dir::String="")
    #restart R session  
    R"""
    sapply(names(sessionInfo()$otherPkgs),function(pkg) detach(paste0('package:',pkg),character.only =T,force = T));rm(list=ls())
    """
    @rput comm_vertices
    @rput ensembl_version
    R"""
    library(biomaRt)
    library(GO.db)
    library(httr)
    library(topGO)
    library(tidyverse)
    
    ## Functional annotation of communities
    vertices = tibble(comm_vertices)
    vertex_names = vertices["label"]
    vertex_names_trimmed = sapply(vertex_names,tools::file_path_sans_ext)
    #get list of string arrays for transcript ids in each community
    comms <- split(vertices,vertices$group)
    comm_ids = lapply(comms,function(x) x$label)
    comm_ids_trimmed = sapply(comm_ids,function(y) sapply(y,function(x) tools::file_path_sans_ext(x)))
    """
    @info "Connecting to biomart ensembl mirror..."
    
    R"""
    ## connect to biomart
    set_config(config(ssl_verifypeer = 0L))
    
    if (ensembl_version=="current")
        {
        ##mirrors to try: "useast" "uswest" "asia"
        ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",mirror="uswest", dataset = "hsapiens_gene_ensembl") 
        go_ids = lapply(comm_ids, function(x) getBM(attributes = c("ensembl_gene_id_version","go_id"),"ensembl_gene_id_version",x,mart = ensembl,useCache = FALSE))
        go_names = vertex_names
        } else
        {
        ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",version=ensembl_version) 
        go_ids = lapply(comm_ids_trimmed, function(x) getBM(attributes = c("ensembl_gene_id","go_id"),"ensembl_gene_id",x,mart = ensembl,useCache = FALSE))
        go_names = vertex_names_trimmed
        }
    """
    @info "Collecting GO terms..."
    R"""
    ##list of all go terms in each community (maybe not necessary with topGO method below now?)
    go_terms = lapply(go_ids,function(x) AnnotationDbi::select(GO.db,columns = "TERM",keys = x$go_id))  
    
    ##Form Gene2GO list required for topGO i.e. gene universe mapped to GO ids
    merged_go_ids = bind_rows(go_ids)
    merged_go_ids[[1]] = factor(merged_go_ids[[1]])
    Genes2GO = split(merged_go_ids[[2]],merged_go_ids[[1]])
    ##list of inputs for topGO, significant genes denoted as those in community (for each community)
    topGOinput = lapply(comm_ids, function (x) factor(as.integer(vertices$label %in% x)))
    #attach IDS
    for (i in 1:length(topGOinput)) {names(topGOinput[[i]]) = go_names} 
    #Generate topGO object for each community
    topGOdata = lapply(topGOinput,function(x) new("topGOdata",ontology = "BP",allGenes = x,annot=annFUN.gene2GO,gene2GO = Genes2GO))
    """
    @info "Running significance tests..."

    R"""
    ##run fisher test on each community
    topGOresults = lapply(topGOdata,function(x) runTest(x,statistic = "fisher"))
    topGOtable = lapply(1:length(topGOdata),function(x) GenTable(topGOdata[[x]],fisher = topGOresults[[x]],ranksOf = "fisher",orderBy="weight"))
    """
    @rget topGOtable
    
    
    if(write_csv==true)
        ##write to CSV
        @info "Writing results to CSV files at $csv_dir..."
        for (key,i) in enumerate(topGOtable)    
            CSV.write("$csv_dir/community_$(key)_annotations.csv",i)
        end
    end
    return topGOtable
end



##bio significance functions
function graphlet_coincidences(vertexlist::Array{String,1},vertex_names::Array{String,1},nametypes::String,adj_matrix::AbstractArray,entrez_id_vector::Array{Int,1},candidates::Dict{String,Array{Int,1}})
        edgelist = NetworkConstruction.edgelist_from_adj(adj_matrix)
        @info "Counting per-edge graphlet relationships..."
        graphlet_counts,Chi,Rel = GraphletCounting.count_graphlets(vertexlist,edgelist,4,run_method="distributed",relationships = true,progress = true)

        ## combine relationships into one array
        rel = vcat(Rel...)
        
        #convert into "nicer" format
        rel_array = broadcast(a->[i for i in a],broadcast(x->x[1:end-1],rel))
        rel_types = last.(rel)
        #clear big rel and Rel TODO do we even need Chi?
        rel = Nothing
        GC.gc()
        ##remove 0 from 3-node entries
        #rel_array = map(y->filter(x->x!=0,y),rel_array)
        #rel_names = map(x->broadcast(y->vertex_gene_names[y],x),rel_array) 
        #rel_transcript_names = map(x->broadcast(y->vertex_names[y],x),rel_array)   
        

        ##split into each homogeneous graphlet and sort there
        ##sort and find unique copies of graphlet 
        @info "Finding unique copies of each graphlet relationship..."
        graphlet_types = string.(unique(last.(split.(collect(keys(graphlet_counts)),"_"))))
        graphlet_rels = Dict{String,Array{Array{Int64,1},1}}()
        @time for g in graphlet_types
            
            #hogs = filter(x->x[end]==g,rel)
            #hogs_array = broadcast(a->[i for i in a],broadcast(x->x[1:end-1],hogs))
            hogs_array = rel_array[findall(x->x==g,rel_types)]

            @info "sorting $g..."
            if(g == "3-path")
                #get rid of leading zero
                hogs_array = map(y->filter(x->x!=0,y),hogs_array)
                graphlet_rels[g] = unique(broadcast(x->[sort(x[[1,3]])[1],x[2],sort(x[[1,3]])[2]],hogs_array))
            end
            if(g == "3-tri")
                #get rid of leading zero
                hogs_array = map(y->filter(x->x!=0,y),hogs_array)
                graphlet_rels[g] = unique(broadcast(x->[sort(x[[1,2,3]])...],hogs_array))
            end
            if(g == "4-path")
                graphlet_rels[g] = unique(broadcast(x->[sort(x[[1,4]])[1],sort(x[[2,3]])...,sort(x[[1,4]])[2]],hogs_array))
            end
            if(g == "4-star")
                graphlet_rels[g] = unique(broadcast(x->[sort(x[[1,2,4]])[1:2]...,x[3],sort(x[[1,2,4]])[3]],hogs_array))
            end
            if(g == "4-tail")
                graphlet_rels[g] = unique(broadcast(x->[sort(x[[1,2]])...,x[3],x[4]],hogs_array))
            end
            if(g == "4-cycle")
                graphlet_rels[g] = unique(broadcast(x->[sort(x[[1,2,3,4]])...],hogs_array))
            end
            if(g == "4-chord")
                graphlet_rels[g] = unique(broadcast(x->[sort(x[[1,4]])[1],sort(x[[2,3]])...,sort(x[[1,4]])[2]],hogs_array))
            end
            if(g == "4-clique")
                graphlet_rels[g] = unique(broadcast(x->[sort(x[[1,2,3,4]])...],hogs_array))
            end
        end
        ## add edges to graphlet_rels as separate entry
        graphlet_rels["2-path"] = [[x...] for x in eachrow(hcat(first.(edgelist),last.(edgelist)))]

#       #graphlet of interest... set manually for now 
#       goi = sig_graphlets.Graphlet[2]     
#       hegoi,hogoi = string.(split(goi,"_")[1:end-1]),string(split(goi,"_")[end])
#       #filter down to homogenous graphlet first
#       hogs = filter(x->x[end]==hogoi,rel)
#       hogs_array = broadcast(a->[i for i in a],broadcast(x->x[1:end-1],hogs))
#       ##sort and find unique copies of graphlet (TODO unique for each graphlet!)
#       if(hogoi == "4-star")
#           hogs_sorted = unique(broadcast(x->[sort(x[[1,2,4]])[1:2]...,x[3],sort(x[[1,2,4]])[3]],hogs_array))
#       end
#   
#       #get those that match heterogeneous pattern as well
#       hogs_types = map(x->broadcast(y->vertexlist[y],x),hogs_sorted)  
#       hegs = hogs_sorted[findall(x->x== hegoi,hogs_types)]
        #get names of transcripts in matching pattern 
#       hegs_names = map(x->broadcast(y->vertex_gene_names[y],x),hegs)  
#       hegs_transcript_names = map(x->broadcast(y->vertex_names[y],x),hegs)    
        #Get KEGG pathway information
        ##now this should be provided as input (to maintain consistency of candidate pathways across functions
        #entrez_id_vector, candidates = get_KEGG_pathways(vertex_names,nametypes) 

        @info "Checking for coincident candidates..."
        Coincidents = Dict{String,Dict{String,Array{Tuple,1}}}()
        for ent in candidates
            @info "checking candidates for $(first(ent))..."
            cands = last(ent)
            #one_coincidents = Array{Array{Int64,1}}(undef,length(keys(graphlet_rels))) 
            #two_coincidents = Array{Array{Int64,1}}(undef,length(keys(graphlet_rels))) 
            #three_coincidents = Array{Array{Int64,1}}(undef,length(keys(graphlet_rels))) 
            #four_coincidents = Array{Array{Int64,1}}(undef,length(keys(graphlet_rels))) 
            #one_coincidents = Array{Tuple,1}()
            two_coincidents = Array{Tuple,1}()
            three_coincidents = Array{Tuple,1}()
            four_coincidents = Array{Tuple,1}()
            for g in keys(graphlet_rels) 
                #graphlets with at least two candidate transcripts involved in process
                #one_coincidents[i] = findall(x->sum(map(y->in(y,x),cands))>0,graphlet_rels[g])
                #two_coincidents[i] = findall(x->sum(map(y->in(y,x),cands))>1,graphlet_rels[g])
                #three_coincidents[i] = findall(x->sum(map(y->in(y,x),cands))>2,graphlet_rels[g])
                #four_coincidents[i] = findall(x->sum(map(y->in(y,x),cands))>3,graphlet_rels[g])
                #push!(one_coincidents,map(x->tuple(graphlet_rels[g][x]...,g),findall(x->sum(map(y->in(y,x),cands))>0,graphlet_rels[g])...))
                push!(two_coincidents,map(x->tuple(graphlet_rels[g][x],g),findall(x->sum(map(y->in(y,x),cands))==2,graphlet_rels[g]))...)
                push!(three_coincidents,map(x->tuple(graphlet_rels[g][x],g),findall(x->sum(map(y->in(y,x),cands))==3,graphlet_rels[g]))...)
                push!(four_coincidents,map(x->tuple(graphlet_rels[g][x],g),findall(x->sum(map(y->in(y,x),cands))==4,graphlet_rels[g]))...)
            end
            #save each in dictionary for candidate
            Coincidents[first(ent)] = Dict("two"=>two_coincidents,"three"=>three_coincidents,"four"=>four_coincidents) 
        end

        #Get data into wide form dataframe, with info on transcript type, ensembl code, entrez id etc...
        Coincidents_df = DataFrame(Pathway=String[],Coincident_type = String[], Hom_graphlet = String[],Vertices = Array{Int64,1}[],Ensembl = Array{String,1}[],Entrez = Array{Int64,1}[],Transcript_type = Array{String,1}[])
        for e in Coincidents
            @info "expanding $(first(e))"
            for ee in last(e)
                for eee in last(ee)
                    push!(Coincidents_df,(first(e),first(ee),last(eee),first(eee),broadcast(x->vertex_names[x],first(eee)),broadcast(x->entrez_id_vector[x],first(eee)),broadcast(x->vertexlist[x],first(eee))))
                end
            end
        end

        ##add inclusion pattern (true if node is in pathway)
        Coincidents_df.Inclusion = [ map(x-> in(x,candidates[Coincidents_df.Pathway[i]]),Coincidents_df.Vertices[i]) for i in 1:size(Coincidents_df)[1]]


        #coincident_names = map(x->broadcast(y->vertex_names[y],x),first.(Coincidents["Morphine addiction"]["three"]))
        #coincident_types = map(x->broadcast(y->vertexlist[y],x),first.(Coincidents["Morphine addiction"]["three"]))
        #coincident_entrez_ids = map(x->broadcast(y->entrez_id_vector[y],x),first.(Coincidents["Morphine addiction"]["three"]))
    return Coincidents_df
end 
function restart_R()

        #restart R session  
        R"""
        sapply(names(sessionInfo()$otherPkgs),function(pkg) detach(paste0('package:',pkg),character.only =T,force = T));rm(list=ls())
        """
end

function biomaRt_connect()
   
    R"""
        library(biomaRt)
        library(httr)
        library(tidyverse)
#       ## connect to biomart
        set_config(config(ssl_verifypeer = 0L))
        ensembl_version = "current" 
        if (ensembl_version=="current")
            {
            ##mirrors to try: "useast" "uswest" "asia"
            ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",mirror="uswest", dataset = "hsapiens_gene_ensembl") 
            } else
            {
            ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",version=ensembl_version) 
            }
        """
end

function get_KEGG_pathways(vertex_names::Array{String,1},nametype::String) 
    restart_R()
    #use small sample for now
    @info "Geting KEGG matches..."
    @rput vertex_names
    biomaRt_connect()
    R"""
        library(edgeR)
        library(tidyverse)
    """
    
    if(nametype == "transcripts")
        R"""
       
            ##Alternative approach: select KEGG terms of interest from whole set of transcripts first, and then look for those in graphlets
            ## full list of entrez_ids mapped to transcripts/genes (neither will be one-to-one, use whichever offers better coverage. Transcripts also preferred as a better reflection of the data rather than expanding to the gene level).
            ##transcripts
            transcripts_trimmed = sapply(vertex_names,tools::file_path_sans_ext)    
            transcripts = data.frame(trimmed_names = transcripts_trimmed)
            entrez_from_transcripts = getBM(attributes=c("ensembl_transcript_id","ensembl_gene_id","entrezgene_id"),"ensembl_transcript_id",transcripts_trimmed, mart = ensembl,useCache=FALSE) 
            transcript_coverage = length(entrez_from_transcripts[[3]])-sum(is.na(entrez_from_transcripts[[3]]))
            transcripts$entrez_id = entrez_from_transcripts[match(transcripts_trimmed,entrez_from_transcripts[[1]]),3]
            
            #get list of entrez ids mapped to KEGG pathways 
            KEGG <-getGeneKEGGLinks(species.KEGG="hsa")
            ## get top hits to select from
            top_terms = topKEGG(kegga(entrez_from_transcripts[[3]],n=Inf,truncate = 34))
            ##get the network candidates for each pathway
            per_pathway = sapply(1:nrow(top_terms),function(x) KEGG$GeneID[ KEGG$PathwayID == row.names(top_terms)[x]])
            in_network = lapply(per_pathway,function(x) transcripts$entrez_id %in% x)
            names(in_network) = top_terms$Pathway
            entrez_id_vector = transcripts$entrez_id
        """ 
    elseif(nametype == "genes")
        R"""
            ##genes
            genes_trimmed = sapply(vertex_names,tools::file_path_sans_ext) 
            genes = data.frame(trimmed_names = genes_trimmed)
            entrez_from_genes = getBM(attributes=c("ensembl_gene_id","entrezgene_id"),"ensembl_gene_id",genes_trimmed, mart = ensembl,useCache=FALSE)   
            gene_coverage = length(entrez_from_genes[[2]])-sum(is.na(entrez_from_genes[[2]]))
            genes$entrez_id = entrez_from_genes[match(genes_trimmed,entrez_from_genes[[1]]),2]
            
            #get list of entrez ids mapped to KEGG pathways 
            KEGG <-getGeneKEGGLinks(species.KEGG="hsa")
            ## get top hits to select from
            top_terms = topKEGG(kegga(entrez_from_genes[[3]],n=Inf,truncate = 34))
            ##get the network candidates for each pathway
            per_pathway = sapply(1:nrow(top_terms),function(x) KEGG$GeneID[ KEGG$PathwayID == row.names(top_terms)[x]])
            in_network = lapply(per_pathway,function(x) genes$entrez_id %in% x)
            names(in_network) = top_terms$Pathway
            entrez_id_vector = genes$entrez_id
        """
    end
    @rget entrez_id_vector
    @rget top_terms
    ##get rid of missing ids (sset to 0 --for now?)
    replace!(entrez_id_vector,missing=>0)
    entrez_id_vector = Int.(entrez_id_vector)
    @rget in_network
    @info "Finding candidates that match top KEGG pathways..."
    candidates = Dict{String,Array{Int,1}}()
    for e in in_network
        candidates[string(first(e))] = findall(.==(true),last(e))
    end
    return (entrez_id_vector, candidates,top_terms)
end

                function pathways_per_node_dict(node_set::Array{Int,1},candidates::Dict{String,Array{Int,1}})
                    output_dict = Dict{Int,Array{String,1}}() 
                    for n in node_set
                        list = []
                        for c in candidates
                            if (n in last(c))
                                push!(list,first(c))
                            end
                        end
                        output_dict[n] = list
                    end
                    return output_dict
                end


function pernode_significance(i::Int,sub_Coincidents::DataFrame,candidate_pathways::Array{String,1},in_key::BitArray{1})#all coincident graphlets that i is involvedF in
    ##for orbit level significance, this function takes a specific node and finds all its coincident graphlets, and then counts the orbit position the node is in in each. returns a dataframe detailing these stats
    ##find those graphlets that include i
    pre_graphlets = filter(:Vertices=> x -> in(i,x),sub_Coincidents)
    #filter further to find cases where at least 2 OTHER nodes in graphlet are in pathway (i.e. self coincidence does not count)
    path_graphlets = filter(:Coincident_type=>x->x!="two",filter(:Pathway=>x-> in(x,candidate_pathways[in_key]),pre_graphlets))
    nonpath_graphlets = filter(:Pathway=>x-> !in(x,candidate_pathways[in_key]),pre_graphlets)
    graphlets = vcat(path_graphlets,nonpath_graphlets) 
    ##set up orbit templates to be checked against
                #peripheral: degree one orbit
                #central: degree two orbit
                #supercentral: degree three orbit
    orbit_templates = Dict("3-path" => Dict(("central" => [0,1,0]), ("peripheral" => [1,0,1])),
                           "3-tri" => Dict(("central" => [1,1,1])),
    "4-path" => Dict(("peripheral" => [1,0,0,1]), ("central" => [0,1,1,0])),
    "4-star" => Dict(("peripheral" => [1,1,0,1]), ("supercentral" => [0,0,1,0])),
    "4-tail" => Dict(("peripheral" => [0,0,0,1]), ("central" => [1,1,0,0]), ("supercentral" => [0,0,1,0])),
    "4-cycle" => Dict(("central"=>[1,1,1,1])),
    "4-chord" => Dict(("supercentral"=>[0,1,1,0]),("central"=>[1,0,0,1])),
    "4-clique" => Dict(("supercentral"=>[1,1,1,1]))) 
     # setup a counter for i for each pathway that features a coincident graphlet of i
     #i_counter = Dict{String,Dict{String,Dict{String,Int64}}}()
    # for p in keys(countmap(graphlets.Pathway))
    #     i_counter[p] =  Dict("3-path" => Dict(("central" => 0), ("peripheral" => 0)),
    # "3-tri" => Dict(("central" => 0)),
    # "4-path" => Dict(("peripheral" => 0), ("central" => 0)),
    # "4-star" => Dict(("peripheral" => 0), ("supercentral" => 0)),
    # "4-tail" => Dict(("peripheral" => 0), ("central" => 0), ("supercentral" => 0)),
    # "4-cycle" => Dict(("central"=>0)),
    # "4-chord" => Dict(("supercentral"=>0),("central"=>0)),
    # "4-clique" => Dict(("supercentral"=>0))) 
    # end
     column =[]
     for g in eachrow(graphlets)
         #find which position i is in this graphlet
         position = g.Vertices.==i
         #find which orbit this position matches (for the given graphlet of g)
         for orb in keys(orbit_templates[g.Hom_graphlet])
             if(Bool(sum(orbit_templates[g.Hom_graphlet][orb].*position)))
                 #i_counter[g.Pathway][g.Hom_graphlet][orb] += 1 
                 push!(column,orb) 
             end
         end
     end
     #append per graphlet label to i's graphlet subset
     graphlets.orbit = column
     # data matrix to tally significance for each pathway in table: (first column peripheral, second column central, third supercentral)
     significance = zeros(Int,length(candidate_pathways),3)
     #per pathway:
     for (i,p) in enumerate(candidate_pathways)
         #collect count of each significance term for this pathway (if the term does not exist, default to 0)
         c = DefaultDict(0,countmap(filter(:Pathway=>x->x==p,graphlets).orbit))
         significance[i,1] = c["peripheral"]
         significance[i,2] = c["central"]
         significance[i,3] = c["supercentral"]
     end
     #pair data with pathway labels into a (per-node) dataframe
     df = DataFrame(Pathway = candidate_pathways, Peripheral = significance[:,1], Central = significance[:,2], Supercentral = significance[:,3])  
     return df
     #@info "Finished $i..."
end

function pernode_significance_detail(i::Int,sub_Coincidents::DataFrame,candidate_pathways::Array{String,1},in_key::BitArray{1})#all coincident graphlets that i is involvedF in
    ##for orbit level significance, this function takes a specific node and finds all its coincident graphlets, and then counts the orbit position the node is in in each. returns a dataframe detailing these stats
    ##find those graphlets that include i
    pre_graphlets = filter(:Vertices=> x -> in(i,x),sub_Coincidents)
    #filter further to find cases where at least 2 OTHER nodes in graphlet are in pathway (i.e. self coincidence does not count)
    path_graphlets = filter(:Coincident_type=>x->x!="two",filter(:Pathway=>x-> in(x,candidate_pathways[in_key]),pre_graphlets))
    nonpath_graphlets = filter(:Pathway=>x-> !in(x,candidate_pathways[in_key]),pre_graphlets)
    graphlets = vcat(path_graphlets,nonpath_graphlets) 
    ##set up orbit templates to be checked against
                #peripheral: degree one orbit
                #central: degree two orbit
                #supercentral: degree three orbit
    orbit_templates = Dict("3-path" => Dict(("central" => [0,1,0]), ("peripheral" => [1,0,1])),
                           "3-tri" => Dict(("central" => [1,1,1])),
    "4-path" => Dict(("peripheral" => [1,0,0,1]), ("central" => [0,1,1,0])),
    "4-star" => Dict(("peripheral" => [1,1,0,1]), ("supercentral" => [0,0,1,0])),
    "4-tail" => Dict(("peripheral" => [0,0,0,1]), ("central" => [1,1,0,0]), ("supercentral" => [0,0,1,0])),
    "4-cycle" => Dict(("central"=>[1,1,1,1])),
    "4-chord" => Dict(("supercentral"=>[0,1,1,0]),("central"=>[1,0,0,1])),
    "4-clique" => Dict(("supercentral"=>[1,1,1,1]))) 
     # setup a counter for i for each pathway that features a coincident graphlet of i
     i_counter = Dict{String,Dict{String,Dict{String,Int64}}}()
     for p in candidate_pathways
         i_counter[p] =  Dict( "4-path" => Dict(("peripheral" => 0), ("central" => 0)),
     "4-star" => Dict(("peripheral" => 0), ("supercentral" => 0)),
     "4-tail" => Dict(("peripheral" => 0), ("central" => 0), ("supercentral" => 0)),
     "4-cycle" => Dict(("central"=>0)),
     "4-chord" => Dict(("supercentral"=>0),("central"=>0)),
     "4-clique" => Dict(("supercentral"=>0))) 
    ##three node version needed
    #"3-path" => Dict(("central" => 0), ("peripheral" => 0)),
    # "3-tri" => Dict(("central" => 0)),
     end

     ##calculate orbit counts across coincident graphlets
     #column =[]
     for g in eachrow(graphlets)
         #find which position i is in this graphlet
         position = g.Vertices.==i
         #find which orbit this position matches (for the given graphlet of g)
         for orb in keys(orbit_templates[g.Hom_graphlet])
             if(Bool(sum(orbit_templates[g.Hom_graphlet][orb].*position)))
                 i_counter[g.Pathway][g.Hom_graphlet][orb] += 1 
                 #push!(column,orb) 
             end
         end
     end
     
     #find score for each individual orbit
     orbit_names = vcat(map(y->map(x->collect(keys(i_counter[candidate_pathways[1]]))[y]*"_"*x, map(x->collect(keys(last(x))),collect(i_counter[candidate_pathways[1]]))[y]),1:6)...)
     orbit_scores = zeros(Int,length(candidate_pathways),length(orbit_names)) 
     for (i,p) in enumerate(candidate_pathways)
         orbit_scores[i,:] = vcat(map(x->collect(values(last(x))),collect(i_counter[p]))...)
     end
     df = DataFrame(orbit_scores,orbit_names)
     insertcols!(df,1,:Pathway=>candidate_pathways)
     return df
     #@info "Finished $i..."
end
