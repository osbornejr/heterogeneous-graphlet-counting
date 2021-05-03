using LightGraphs,RCall
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
	x = length(unique(communities$membership))-11
	colour_palette = c(brewer.pal(name = "Spectral", n = 11),brewer.pal(name = "BrBG", n = x))
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




function get_functional_annotations(comm_vertices::DataFrame;ensembl_version::String="current",write_csv::Bool = true,csv_dir::String)
	
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
