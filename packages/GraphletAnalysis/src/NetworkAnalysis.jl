using Graphs,RCall,DataFrames,GraphletCounting,NetworkConstruction
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


function get_community_structure(adj_matrix::AbstractArray,vertex_names::Vector{<:AbstractString},detection_method::String;plot_prefix::String = "", threejs_plot::Bool = false)
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
        ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",mirror="useast", dataset = "hsapiens_gene_ensembl") 
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
function per_graphlet_coincidence(graphlet::AbstractArray,type::AbstractString,candidates)

end
function graphlet_coincidences(rel::Matrix{Int},rel_types::AbstractVector,vertexlist::Vector{<:AbstractString},vertex_names::Vector{<:AbstractString},entrez_id_vector::Array{Int,1},candidates::Dict{String,Array{Int,1}})
        #convert into "nicer" format
        #rel_array = eachrow(rel)
        
        #rel = Nothing
        ##remove 0 from 3-node entries
        #rel_array = map(y->filter(x->x!=0,y),rel_array)
        #rel_names = map(x->broadcast(y->vertex_gene_names[y],x),rel_array) 
        #rel_transcript_names = map(x->broadcast(y->vertex_names[y],x),rel_array)   
        


        @info "Checking for coincident candidates..."
        Coincidents = Dict{String,Dict{Int,Array{Tuple,1}}}()
        for ent in candidates
            @info "checking candidates for $(first(ent))..."
            cands = last(ent)
            #one_coincidents = Array{Array{Int64,1}}(undef,length(keys(graphlet_rels))) 
            #two_coincidents = Array{Array{Int64,1}}(undef,length(keys(graphlet_rels))) 
            #three_coincidents = Array{Array{Int64,1}}(undef,length(keys(graphlet_rels))) 
            #four_coincidents = Array{Array{Int64,1}}(undef,length(keys(graphlet_rels))) 
            one_coincidents = Array{Tuple,1}()
            two_coincidents = Array{Tuple,1}()
            three_coincidents = Array{Tuple,1}()
            four_coincidents = Array{Tuple,1}()

            #per graphlet
            
            for (i,r) in enumerate(eachrow(rel))

                if(sum(map(x->x in cands,r))==1)         
                    push!(one_coincidents,tuple(r,rel_types[i]))                
                elseif(sum(map(x->x in cands,r))==2)         
                    push!(two_coincidents,tuple(r,rel_types[i]))                
                
                elseif(sum(map(x->x in cands,r))==3)         
                    push!(three_coincidents,tuple(r,rel_types[i]))                

                elseif(sum(map(x->x in cands,r))==4)         
                    push!(four_coincidents,tuple(r,rel_types[i]))                
                end

            end
           
            


            #for g in keys(graphlet_rels) 
            #    #graphlets with at least two candidate transcripts involved in process
            #    #one_coincidents[i] = findall(x->sum(map(y->in(y,x),cands))>0,graphlet_rels[g])
            #    #two_coincidents[i] = findall(x->sum(map(y->in(y,x),cands))>1,graphlet_rels[g])
            #    #three_coincidents[i] = findall(x->sum(map(y->in(y,x),cands))>2,graphlet_rels[g])
            #    #four_coincidents[i] = findall(x->sum(map(y->in(y,x),cands))>3,graphlet_rels[g])
            #    #push!(one_coincidents,map(x->tuple(graphlet_rels[g][x]...,g),findall(x->sum(map(y->in(y,x),cands))>0,graphlet_rels[g])...))
            #    push!(two_coincidents,map(x->tuple(graphlet_rels[g][x],g),findall(x->sum(map(y->in(y,x),cands))==2,graphlet_rels[g]))...)
            #    push!(three_coincidents,map(x->tuple(graphlet_rels[g][x],g),findall(x->sum(map(y->in(y,x),cands))==3,graphlet_rels[g]))...)
            #    push!(four_coincidents,map(x->tuple(graphlet_rels[g][x],g),findall(x->sum(map(y->in(y,x),cands))==4,graphlet_rels[g]))...)
            #end
            #save each in dictionary for candidate
            Coincidents[first(ent)] = Dict(1=>one_coincidents,2=>two_coincidents,3=>three_coincidents,4=>four_coincidents) 
        end

        #Get data into wide form dataframe, with info on transcript type, ensembl code, entrez id etc...
        Coincidents_df = DataFrame(Pathway=String[],Coincident_nodes = Int[], Hom_graphlet = String[],Vertices = Array{Int64,1}[],Ensembl = Array{String,1}[],Entrez = Array{Int64,1}[],Transcript_type = Array{String,1}[])
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
#       ## connect to biomart
        set_config(config(ssl_verifypeer = 0L))
        ensembl_version = "current" 
        if (ensembl_version=="current")
            {
            ##mirrors to try: "useast" "uswest" "asia"
            ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",mirror="useast", dataset = "hsapiens_gene_ensembl") 
            } else
            {
            ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",version=ensembl_version) 
            }
        """
end

"""
    get_entrez_ids(names,nametype)

Generate a map between `names` and the entrez gene ids associated with them. 
Elements of `names` may be any string, but only ensembl names will generate a match. 
"Name.x" will have ".x" automatically trimmed as biomaRt will not match to the sub id level.

Matches be either at the "transcript" or "gene" level, which is given by `nametype`.

"""
function get_entrez_ids(names::Vector{<:AbstractString},nametype::String)

    restart_R()
    biomaRt_connect()
    @info "Getting Entrez ids..."
    if (nametype == "transcript")
        ensembl_search_term = "ensembl_transcript_id"
    elseif (nametype=="gene")
        ensembl_search_term = "ensembl_gene_id"
    else
        throw(ArgumentError("nametype must be either 'transcript' or 'gene'."))
    end

    @rput names
    @rput ensembl_search_term 
    R"""

    ## full list of entrez_ids mapped to names (either transcript or gene names, neither will be one-to-one, use whichever offers better coverage. Transcripts also preferred as a better reflection of the data rather than expanding to the gene level).
    
    ### we trim names to remove transcript specific reference (better for ensembl hits)
    names_trimmed = sapply(names,tools::file_path_sans_ext)    
    ##create dataframe to store map between original names and entrez ids
    #(form dataframe with trimmed names)
    name_map = data.frame(name =names,trimmed_name = names_trimmed)
    
    #this will map each name to an entrez id, IF the name exists in ensembl database 
    entrez_from_names = getBM(attributes=c("ensembl_transcript_id","ensembl_gene_id","entrezgene_id"),ensembl_search_term,names_trimmed, mart = ensembl,useCache=FALSE) 
    ##add column denoting if name was found in ensembl database
    name_map$ensembl_match = !is.na(match(names_trimmed,entrez_from_names[[ensembl_search_term]]))
    ## add FIRST entrez id match to map dataframe 
    name_map$entrez_id = entrez_from_names[match(names_trimmed,entrez_from_names[[ensembl_search_term]]),3]
    name_coverage = length(entrez_from_names[[3]])-sum(is.na(entrez_from_names[[3]]))
    """
    @rget name_map
    ##append full set of associated entrez ids (name_map will only have FIRST match
    @rget entrez_from_names
    entrez_set(x) = unique(filter(ensembl_search_term => ==(x),entrez_from_names).entrezgene_id)
    name_map.entrez_ids = map(x->entrez_set(x),name_map.trimmed_name)
    return name_map
end

function get_GO_terms(vertex_names::Vector{<:AbstractString},nametype::String)

    ##get name map to entrez ids for vertex names
    ## note: R restart and biomaRt init happen in this function, so dont need to do again here
    name_map = get_entrez_ids(vertex_names,nametype)
    
    ##just select one entrez id (first match) atm, and remove transcripts with no corresponding entrez id
    entrez_map = filter(:entrez_id=>!ismissing,name_map)
    ## importantly, record which nodes in network ahve carried through
    entrez_map.network_node_id = findall(!ismissing,name_map.entrez_id)
    entrez_ids = convert(Vector{Int},entrez_map.entrez_id)
    #restart_R()
    @info "Getting GO matches..."
    @rput entrez_ids
    #biomaRt_connect()
    R"""
    library(edgeR)

    ## get top hits to select from
    top_terms = topGO(goana(entrez_ids))
    """
    @rget top_terms 

    return top_terms
end

function get_KEGG_pathways(vertex_names::Vector{<:AbstractString},nametype::String) 
    
    ##get name map to entrez ids for vertex names
    ## note: R restart and biomaRt init happen in this function, so dont need to do again here
    name_map = get_entrez_ids(vertex_names,nametype)
    
    ##just select one entrez id (first match) atm, and remove transcripts with no corresponding entrez id
    entrez_map = filter(:entrez_id=>!ismissing,name_map)
    ## importantly, record which nodes in network ahve carried through
    entrez_map.network_node_id = findall(!ismissing,name_map.entrez_id)
    entrez_ids = convert(Vector{Int},entrez_map.entrez_id)
    #restart_R()
    @info "Getting KEGG matches..."
    @rput entrez_ids
    #biomaRt_connect()
    R"""
    library(edgeR)

    #get list of entrez ids mapped to KEGG pathways 
    ##For some reason this is not working atm (SSH issue? MAybe Kitty? TODO) so we will manually load in pathway info
    #KEGG <-getGeneKEGGLinks(species.KEGG="hsa")
    pathway_links <- read.table("data/kegg_pathway_links_hsa.txt")
    names(pathway_links) <- c("GeneID","PathwayID")
    ##need to convert gene ids to correct integer format
    pathway_links$GeneID = strtoi(sapply(pathway_links$GeneID,function(x) sub("hsa:","",x)))
    pathway_list <- read.table("data/kegg_pathway_list_hsa.txt",sep = "\t")
    names(pathway_list) <- c("PathwayID","PathwayName")

    ## get top hits to select from
    top_terms = topKEGG(kegga(entrez_ids,n=Inf,truncate = 34,gene.pathway = pathway_links,pathway.names = pathway_list))

    """ 
    @rget top_terms
    return top_terms
end

function get_KEGG_candidates(vertex_names::Vector{<:AbstractString},nametype::String)
    ##get name map to entrez ids for vertex names
    ## note: R restart and biomaRt init happen in this function, so dont need to do again here
    name_map = get_entrez_ids(vertex_names,nametype)
    
    ##just select one entrez id (first match) atm, and remove transcripts with no corresponding entrez id
    entrez_map = filter(:entrez_id=>!ismissing,name_map)
    ## importantly, record which nodes in network ahve carried through
    entrez_map.network_node_id = findall(!ismissing,name_map.entrez_id)
    entrez_ids = convert(Vector{Int},entrez_map.entrez_id)
    #run this to get top_terms
    top_terms = get_KEGG_pathways(vertex_names,nametype)

    @rput entrez_ids
    @rput top_terms 
    R"""
    #get list of entrez ids mapped to KEGG pathways 
    ##For some reason this is not working atm (SSH issue? MAybe Kitty? TODO) so we will manually load in pathway info
    #KEGG <-getGeneKEGGLinks(species.KEGG="hsa")
    pathway_links <- read.table("data/kegg_pathway_links_hsa.txt")
    names(pathway_links) <- c("GeneID","PathwayID")
    ##need to convert gene ids to correct integer format
    pathway_links$GeneID = strtoi(sapply(pathway_links$GeneID,function(x) sub("hsa:","",x)))
    pathway_list <- read.table("data/kegg_pathway_list_hsa.txt",sep = "\t")
    names(pathway_list) <- c("PathwayID","PathwayName")
    
    ##get the network candidates for each pathway
    per_pathway = sapply(1:nrow(top_terms),function(x) pathway_links$GeneID[ pathway_links$PathwayID == row.names(top_terms)[x]])
    in_network = lapply(per_pathway,function(x) entrez_ids %in% x)
    names(in_network) = top_terms$Pathway
    """
    @rget in_network
    @info "Finding candidates that match top KEGG pathways..."
    candidates = Dict{String,Array{Int,1}}()
    for e in in_network
        candidates[string(first(e))] = findall(.==(true),last(e))
    end
    #TODO check that reduction to only those transcripts with a entrez match here does not cause problems down the line (coincident analysis) 
    return (entrez_map,candidates,top_terms)
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
    ### NOTE: Currently depreceated. Will need to be updated in line with detail function/ 
    ##for orbit level significance, this function takes a specific node and finds all its coincident graphlets, and then counts the orbit position the node is in in each. returns a dataframe detailing these stats
    ##find those graphlets that include i
    pre_graphlets = filter(:Vertices=> x -> in(i,x),sub_Coincidents)
    #filter further to find cases where at least 2 OTHER nodes in graphlet are in pathway (i.e. self coincidence does not count)
    path_graphlets = filter(:Coincident_nodes=>x->x>2,filter(:Pathway=>x-> in(x,candidate_pathways[in_key]),pre_graphlets))
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

function pernode_significance_detail(i::Int,sub_Coincidents::DataFrame,graphlet_size::Int,candidate_pathways::Array{String,1},in_key::BitArray{1})#all coincident graphlets that i is involvedF in
    ##for orbit level significance, this function takes a specific node and finds all its coincident graphlets, and then counts the orbit position the node is in in each. returns a dataframe detailing these stats
    ##find those graphlets that include i
    pre_graphlets = filter(:Vertices=> x -> in(i,x),sub_Coincidents)
    #filter further to find cases where at least `thresh` OTHER nodes in graphlet are in pathway (i.e. self coincidence does not count)
    #TODO coincidence threshold settable in config? May need logic check against graphlet size
    thresh = graphlet_size - 1
    path_graphlets = filter(:Coincident_nodes=>x->x>thresh,filter(:Pathway=>x-> in(x,candidate_pathways[in_key]),pre_graphlets))
    ##for nonpath, we only require coincident node count to be equal to thresh
    nonpath_graphlets = filter(:Coincident_nodes=>x->x>=thresh,filter(:Pathway=>x-> in(x,candidate_pathways[in_key]),pre_graphlets))
    graphlets = vcat(path_graphlets,nonpath_graphlets) 
    ##set up orbit templates to be checked against
                #peripheral: degree one orbit
                #central: degree two orbit
                #supercentral: degree three orbit
                ##all possible 3 and 4 nodes)
                ## note: we set up as a template here to ensure that all orbits are recorded, not just those that occur in node (or even network) subset 
    if (graphlet_size == 3)
        orbit_templates = Dict("3-path" => Dict(("central" => [0,1,0]), ("peripheral" => [1,0,1])),
                           "3-tri" => Dict(("central" => [1,1,1])))
        # setup a counter for i for each pathway that features a coincident graphlet of i
        empty_counter = Dict("3-path" => Dict(("central" => 0), ("peripheral" => 0)),
                          "3-tri" => Dict(("central" => 0)))
    elseif (graphlet_size == 4)
    
        orbit_templates = Dict("4-path" => Dict(("peripheral" => [1,0,0,1]), ("central" => [0,1,1,0])),
                           "4-star" => Dict(("peripheral" => [1,1,0,1]), ("supercentral" => [0,0,1,0])),
                           "4-tail" => Dict(("peripheral" => [0,0,0,1]), ("central" => [1,1,0,0]), ("supercentral" => [0,0,1,0])),
                           "4-cycle" => Dict(("central"=>[1,1,1,1])),
                           "4-chord" => Dict(("supercentral"=>[0,1,1,0]),("central"=>[1,0,0,1])),
                           "4-clique" => Dict(("supercentral"=>[1,1,1,1]))) 
     # setup a counter for i for each pathway that features a coincident graphlet of i
     empty_counter = Dict("4-path" => Dict(("peripheral" => 0), ("central" => 0)),
                           "4-star" => Dict(("peripheral" => 0), ("supercentral" => 0)),
                           "4-tail" => Dict(("peripheral" => 0), ("central" => 0), ("supercentral" => 0)),
                           "4-cycle" => Dict(("central"=>0)),
                           "4-chord" => Dict(("supercentral"=>0),("central"=>0)),
                           "4-clique" => Dict(("supercentral"=>0))) 
    
    end
    ## store all pathway counters in this dict
     i_counter = Dict{String,Dict{String,Dict{String,Int64}}}()
     for p in candidate_pathways 
         i_counter[p] =  deepcopy(empty_counter)
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
     
     #find score for each individual orbit TODO fix reliance on 1:8 vector covering arbitrary number of graphlet entries
     orbit_names = vcat(map(y->map(x->collect(keys(i_counter[candidate_pathways[1]]))[y]*"_"*x, map(x->collect(keys(last(x))),collect(i_counter[candidate_pathways[1]]))[y]),1:length(orbit_templates))...)
     orbit_scores = zeros(Int,length(candidate_pathways),length(orbit_names)) 
     for (i,p) in enumerate(candidate_pathways)
         orbit_scores[i,:] = vcat(map(x->collect(values(last(x))),collect(i_counter[p]))...)
     end
     df = DataFrame(orbit_scores,orbit_names)
     insertcols!(df,1,:Pathway=>candidate_pathways)
     return df
     #@info "Finished $i..."
end

##ANalysis of coincidents (messy and incoherent) TODO needs to be formalised and put into proper package places
    ##find which types are excluded in general, and then only as a cause of having no Entrez id
    #Coincidents.excluded = [Coincidents.Transcript_type[i][Coincidents.Inclusion[i].==0] for i in 1:size(Coincidents)[1]]
    #Coincidents.excluded_Entrez = [Coincidents.Transcript_type[i][Coincidents.Entrez[i].==0] for i in 1:size(Coincidents)[1]]
    ##find only those coincidents that involve non-coding transcripts
    #Coincidents_noncoding = Coincidents[findall(x-> "noncoding" in x, Coincidents.Transcript_type),:]

    ##Nonuniformity test: finds graphlets where the included nodes differ across different pathways (sampling for now for speed)
    #graphlet = "4-star"
    #nonuniforms = []
    #for y in filter(p->last(p)>1,countmap(filter(:Hom_graphlet=> x-> x == graphlet, Coincidents[1:199000,:]).Vertices))
    #    test = sum(filter(:Vertices => x-> x == first(y),filter(:Hom_graphlet => x-> x == graphlet,Coincidents))[!,8])/last(y)
    #    if (sum(((test.>0) - (test.<1)).==0)>0)
    #        push!(nonuniforms,first(y))
    #    end
    #end
    ##types of exlusions: which transcript types are most likely to be missing from the pathway in a graphlet
    #countmap([Coincidents.Transcript_type[i][Coincidents.Inclusion[i].==0] for i in 1:size(Coincidents)[1]])
    #countmap([Coincidents_noncoding.Transcript_type[i][Coincidents_noncoding.Inclusion[i].==0] for i in 1:size(Coincidents_noncoding)[1]])

    ##agreement between entrez and inclusion info
    #sum(map(x->x.!==0,Coincidents.Entrez).==Coincidents.Inclusion)
    #sum(map(x->x.!==0,Coincidents_noncoding.Entrez).==Coincidents_noncoding.Inclusion)


 

   # # this table stores the average for each pathway (of nodes that are known to be in the pathway)
  #  significance_bars = zeros(length(keys(candidates)),size(orbit_sigs_array[1])[2])
  #  for (i,c) in enumerate(keys(candidates))
  #      significance_bars[i,:] = (sum(map(x->Array(x[!,2:end]),orbit_sigs[candidates[c]]))./length(vertexlist))[i,:]
  #  end
  #  ##now compare bars against profiles of non-pathway nodes
  #  ## choose a subset of nodes to look at. Can be boolean BitArray (with length equal to all nodes) or a specific list of nodes 
  #  subset = entrez_id_vector.==0
  #  putative_pathways = Array{Array{String,1}}(undef,length(orbit_sigs_array[subset]))
  #  for (i,t) in enumerate(orbit_sigs_array[subset])
  #      putative_pathways[i] = collect(keys(candidates))[vec((sum(t.>significance_bars,dims=2).>2))]
  #  end

    # Gadfly beeswarm visualisation:
 #   output_dir = "$cwd/output/plots/orbit_significance_$(orbit_sigs_method)/"
 #   run(`mkdir -p $output_dir`)
 #   # get data into wide format
 #   # size of each df
 #   last_col = size(orbit_sigs[1])[2]-1
 #   wide_orbit_sigs = vcat(map(x->stack(x,2:last_col+1),orbit_sigs)...)
 #   #for each pathway, we map the three orbit categories side by side
 #   palette = ["#db63c5","#bababa","#32a852"]
#    for p in keys(candidates)
#        @info "Drawing beeswarm for $p..."    
#        p_df = filter(:Pathway=>x->x == p,wide_orbit_sigs)
#        #insertcols!(p_df,:log_value =>log.(p_df.value))
#        #define colors by whether entrez id of node is in pathway
#        in_pathway = vcat(collect(eachrow(repeat(in.(1:length(vertexlist),Ref(candidates[p])),1,last_col)))...)
#        non_entrez = vcat(collect(eachrow(repeat(entrez_id_vector.==0,1,last_col)))...)
#        coloring = CategoricalArray(in_pathway-non_entrez)
#        coloring = recode(coloring,-1=>"unidentified",0=>"not in pathway",1=>"in pathway")
#        insertcols!(p_df,:color =>coloring)
#        #remove non pathway nodes
#        #filter!(:color=>x->x!="not in pathway",p_df)
#        #remove zero nodes
#        filter!(:value=>x->x!=0,p_df)
#        pl = plot(p_df,x = :variable,y = :value, color = :color,Guide.title(p),Geom.beeswarm(padding = 1mm),Theme(bar_spacing=1mm,point_size=0.5mm),Scale.color_discrete_manual(palette...));
#        draw(SVG("$(output_dir)/$(p)_beeswarm.svg",30cm,20cm),pl)
#    end
    #@time motif_counts = find_motifs(edgelist,"hetero_rewire",100, typed = true, typelist = vec(vertexlist),plotfile="$cache_dir/motif_detection.svg",graphlet_size = 4)

    #High zero exploration
    #TODO remove this and just use low_filter method below to more organically acheive same thing
    ## over all transcripts 
   # total_zero_proportion = sum(map(x->x.==0,orbit_sigs_array))./length(orbit_sigs_array)
   # ##over a subset
   # subset = candidates[candidate_pathways[1]] 
   # subset_zero_proportion = sum(map(x->x.==0,orbit_sigs_array[subset]))./length(orbit_sigs_array[subset])
   # #compare
   # zero_comparison = subset_zero_proportion.<total_zero_proportion
   # #compare for all pathways,speficically for that pathway
   # zero_scores = zeros(Int,length(candidate_pathways),last_col) 
   # for (i,p) in enumerate(candidate_pathways)
   #     subset = candidates[p] 
   #     subset_zero_proportion = sum(map(x->x.==0,orbit_sigs_array[subset]))./length(orbit_sigs_array[subset])
   #     #compare
   #     zero_scores[i,:] = (subset_zero_proportion.<total_zero_proportion)[i,:]
   # end
   # #find those pathways with majority of zero proportions below total proportions 
   # zero_passes = vec(sum(zero_scores,dims=2).>(last_col/2))
   # zero_candidate_pathways = candidate_pathways[zero_passes]
   # zero_orbit_sigs = map(x->filter(:Pathway=>y->y in zero_candidate_pathways,x),orbit_sigs)
   # zero_orbit_sigs_array = map(x->Array(x[!,2:end]),zero_orbit_sigs)
   # zero_candidates = Dict(Pair.(zero_candidate_pathways,[candidates[x] for x in zero_candidate_pathways]))

   # #uniqueness of pathway contributors (concerns of too much overlap) 
   # sig_pathway_occurences = countmap(vcat([candidates[x] for x in zero_candidate_pathways]...))
   # m = max(collect(values(sig_pathway_occurences))...) 
   # #m = 8
   # supersharers = first.(filter(x->last(x)==m,collect(sig_pathway_occurences)))
   # #for these supersharers, find the set of pathways they are involved in
   # supersharer_pathways = GraphletAnalysis.pathways_per_node_dict(supersharers,zero_candidates)
   # in_group = collect(keys(countmap(vcat(collect(values(supersharer_pathways))...))))
   # not_in_group = zero_candidate_pathways[.!(in.(zero_candidate_pathways,Ref(collect(keys(countmap(vcat(collect(values(supersharer_pathways))...)))))))]
   # countmap(collect(values(supersharer_pathways)))
   # #do we need to rule out pathways dominated by supersharers? TODO


#    #ecdfs
#    #store ecdf functions in table
#    ecdf_table = Array{ECDF,2}(undef,length(zero_candidate_pathways),last_col)
#    for (i,p) in enumerate(zero_candidate_pathways)
#        for j in 1:last_col
#            ecdf_table[i,j] = ecdf(map(x->x[i,j],zero_orbit_sigs_array))
#        end
#    end
#    #first check candidate probabilities across all categories
#    known_pathway_probs = Array{Array{Float64,2},1}(undef,length(zero_candidate_pathways))
#    for (i,c) in enumerate(zero_candidate_pathways)
#        known_pathway_probs[i] = hcat([map(ecdf_table[i,j],map(x->x[i,j],orbit_sigs_array[zero_candidates[c]])) for j in 1:last_col]...)
#    end
#
#    #define orbit names here for now
#    orbit_names = names(orbit_sigs[1])[2:last_col+1]
#
#    #collect known pathway vectors for corresponding known pathway nodes
#    known_pathway_dfs = Array{DataFrame,1}(undef,length(zero_candidate_pathways))
#    for (i,p) in enumerate(zero_candidate_pathways)
#        subset = zero_candidates[p]
#        df_build = DataFrame(hcat(map(x->x[i,:],zero_orbit_sigs_array[subset])...)',orbit_names)
#        insertcols!(df_build,1,:shared=>[sig_pathway_occurences[x] for x in subset].-1)
#        insertcols!(df_build,1,:transcript_id=>subset)
#        known_pathway_dfs[i] = df_build 
#        #print("Known nodes in pathway $p are in this many other pathways:\n")
#        #print([sig_pathway_occurences[x] for x in subset].-1)
#        #print("\n")
#    end
#    known_pathway_arrays = map(x->Array(x[!,3:end]),known_pathway_dfs)
#    #known ecdfs
#    #calculate ecdfs just for known pathway node values (each row corresponds to a pathway, each column to an orbit) 
#    known_ecdf_table = Array{ECDF,2}(undef,length(zero_candidate_pathways),last_col)
#    for (i,p) in enumerate(zero_candidate_pathways)
#        known_ecdf_table[i,:] = ecdf.(eachcol(known_pathway_arrays[i]))
#    end
#
#    ##compare each unknown node to known ecdfs
#    #set threshold for which to check on i.e. is node orbit count higher than the top ~thresh~% of known nodes
#    thresh = 0.05
#    ub = 1 -thresh
#    #first check if there are any orbits for any pathways that have such low representation that 0 count exceed threshold. we will filter these out of analysis
#    low_filter = .!(map(x->x(0),known_ecdf_table).>ub)
#    #for each node, map each known ecdf to the corresponding orbit count and check against threshold, factoring in the low filter
#    unknown_ecdf_comparison = map(y->(reshape(map(x->known_ecdf_table[x](y[x]),1:length(known_ecdf_table)),size(known_ecdf_table)[1],size(known_ecdf_table)[2]).>ub).*low_filter,zero_orbit_sigs_array)
#    #determine a node as significantly linked to a pathway if at least half its orbit counts are significant
#    sig_check = map(x->(sum(x,dims=2).>(last_col/2)),unknown_ecdf_comparison)
#    sig_nodes= findall(x->sum(x)>0,sig_check)
#    sig_nodes_dict = Dict(Pair.(sig_nodes,map(x->zero_candidate_pathways[vec(x)],sig_check[sig_nodes])))
#
#    #shape plots into a grid
#    ncols = 6
#    dims = fldmod(length(zero_candidate_pathways),ncols)
#    plots = Array{Union{Plot,Context},2}(undef,dims[1]+(dims[2]>0),ncols)
#    #plots = Array{Union{Plot,Context},1}(undef,length(keys(candidates)))
#    for (i,p) in enumerate(zero_candidate_pathways)
#        #For total ecdfs:
#        #find max over all measures for pathway
#        #                    m = max([map(x->x[i,1],orbit_sigs_array)...,map(x->x[i,2],orbit_sigs_array)..., map(x->x[i,3],orbit_sigs_array)...]...)
#        #                    plots[i] = plot([layer(x->ecdf(map(x->x[i,j],orbit_sigs_array))(x),0,m,color=[j]) for j in 1:last_col]...,
#        #                              Scale.color_discrete_manual("orange", "green", "purple"),
#        #                              Guide.title(p),
#        #                              Guide.xlabel("count"),
#        #                              Theme(major_label_font_size=4pt,key_position=:none));
#        #                              #Guide.colorkey(title="orbit position"),
#        #                              #Guide.title(p));
#        #                              #
#        #for known ecdfs only:
#        m = max(known_pathway_arrays[i]...)
#        plots[i] = plot([layer(x->known_ecdf_table[i,j](x),0,m,color=[orbit_names[j]]) for j in 1:last_col]...,
#                        Scale.color_discrete_manual("orange", "green", "purple"),
#                        Guide.title(p),
#                        Guide.xlabel("count"),
#                        Theme(major_label_font_size=4pt,key_position=:none));
#        #Guide.colorkey(title="orbit position"),
#        #Guide.title(p));
#    end
#    #Add legend pane
#    legend = plot(wide_orbit_sigs,color=:variable,
#                  Geom.blank,
#                  Scale.color_discrete_manual("orange", "green", "purple"));
#    #append blank gridspots if necessary
#    for i in 1:length(plots)
#        if(!isassigned(plots,i))
#            plots[i] = context()
#        end
#    end
#    #TODO this needs to be more generalised for any dimension of gridstack/number of pathways 
#    plots[2,6] = legend;
#    draw(SVG("$(output_dir)/known_ecdfs.svg",30cm,20cm),gridstack(plots))
#
#
#    ## make a table to analyse different candidate levels for coincident graphlets
#
#
#    ##function to get the n permutations of a set xs
#    all_perm(xs, n) = vec(map(collect, Iterators.product(ntuple(_ -> xs, n)...)))
#    ## method to find all orbit permutations of a type set
#    set = ["coding","noncoding"]
#
#    # specify orbit classes
#    orbit = [1,1,1,1]
#    order = length(orbit)
#    ##get all possible combinations from set (before sorting)
#    base = collect.(vcat(collect(Iterators.product(repeat([set],order)...))...))
#    #sort via each orbit equivalence class 
#    for o in unique(orbit)
#        template = orbit.==o
#        for i in base
#            i[template] = sort(i[template])
#        end
#    end
#    #find unique permutations
#    gs = unique(base)  
#
    #return sig_nodes_dict


