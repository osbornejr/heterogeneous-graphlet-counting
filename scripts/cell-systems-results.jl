using ProjectFunctions,CSV,DataStructures,DataFrames


function tex_typed_results_bar(data::DataFrame,outfile::String)
        n_graphlet_types = size(data)[1]
        n_graphlets = sum(data.Observed)
        tex = "\\begin{tikzpicture}
        \\begin{axis}[
        ybar=0pt,% space of 0pt between adjacent bars
        bar width=2,
        width=7cm,
        height=12cm,
        minor x tick num=4,
        xtick=data,
        enlargelimits=0.15,
        ]
        "
        ##add coordinates for bars
        tex*= "\\addplot coordinates {"
        for i in 1:n_graphlet_types
            tex*="($(i),$(data[i,2]/n_graphlets)) "
        end
        tex*="};\n"
        tex*= "\\addplot coordinates {"
        for i in 1:n_graphlet_types
            tex*="($(i),$(data[i,3]/n_graphlets)) "
        end
        tex*="};\n\\end{axis}\n\\end{tikzpicture}"
        write(outfile,tex)
end


function tex_enrichment_bar(top_names::Vector{String},enrichment_scores::Vector{Float64},outfile::String)

        tex = "\\begin{tikzpicture}
        \\begin{axis}[
            xbar,
            xmin = 0,
            width = 5cm,
            height = 10cm,
            enlarge y limits=0.05,
            axis lines* = left,
            xlabel = Enrichment score \$(-\\log(p))\$,
            symbolic y coords ={"
##add in names for y coords
#tex*= join(sort(kegg_top_terms,:enrichment_score).Pathway,",")*"},\n"
tex*= join(reverse(top_names),",")*"},\n"
tex*="ytick=data]\n"
##add coordinates for bars
tex*= "\\addplot coordinates {"
for i in 1:length(top_names)
    tex*="($(enrichment_scores[i]),$(top_names[i])) "
end
tex*="};\n\\end{axis}\n\\end{tikzpicture}"
write(outfile,tex)
end

function tex_degree_distribution(data::Vector{Int},outfile::String;title::String="")
    tex= "\\begin{tikzpicture}
    \\begin{axis}[title=$(title),xlabel=Degree,ylabel=Frequency]
            \\addplot+ [
                       hist={
                             bins=50,
                             handler/.style={sharp plot},
                             intervals=false,
                            },
                            mark=none,
                      ] table [row sep=\\\\,y index=0] {
                                                      data\\\\"
                                                      tex*=join(data,"\\\\")*"\\\\" 
    tex*=                                                 "};
        \\end{axis}
    \\end{tikzpicture}"
    write(outfile,tex)
end


## RESULTS TABLES



#NETWORK STATS
## load in info
raw_counts,processed_counts,similarity_matrix,adj_matrix,network_counts,vertexlist,edgelist = ProjectFunctions.get_output_data()

### Degree distributions
dd = GraphletCounting.typed_degree_distribution(vertexlist,edgelist)
plot_dd_data = DataFrame(coding = map(x->DefaultDict(0,x)["coding"],dd),noncoding = map(x->DefaultDict(0,x)["noncoding"],dd));

##generate degree distribution plots for each edge type
tex_degree_distribution(plot_dd_data[vertexlist.=="coding",:].coding,"output/share/degree-dist-coding-coding.tex",title="Coding degree distribution for coding nodes")
tex_degree_distribution(plot_dd_data[vertexlist.=="coding",:].noncoding,"output/share/degree-dist-coding-noncoding.tex",title="Noncoding degree distribution for coding nodes")
tex_degree_distribution(plot_dd_data[vertexlist.=="noncoding",:].coding,"output/share/degree-dist-noncoding-coding.tex",title="Coding degree distribution for noncoding nodes")
tex_degree_distribution(plot_dd_data[vertexlist.=="noncoding",:].noncoding,"output/share/degree-dist-noncoding-noncoding.tex",title="Noncoding degree distribution for noncoding nodes")

#BIOLOGICAL VALIDATION
##load in information
kegg_top_terms,go_top_terms = ProjectFunctions.biological_validation()
### KEGG enrichment
##bio validation table of top kegg terms associated with network
## calculate enrichment score (invert p-values on log scale)
kegg_top_terms.enrichment_score = broadcast(x->-log(x),kegg_top_terms.P_DE) 
## remove homo sapiens tail from every pathway (mention in caption)
kegg_top_terms.Pathway = replace.(kegg_top_terms.Pathway," - Homo sapiens (human)"=>"")
#rekmove other longwinded subcats
kegg_top_terms.Pathway = replace.(kegg_top_terms.Pathway," - multiple diseases"=>"")
kegg_top_terms.Pathway = replace.(kegg_top_terms.Pathway," - reactive oxygen species"=>"")
kegg_top_terms.Pathway = replace.(kegg_top_terms.Pathway," and "=>"/")

tex_enrichment_bar(kegg_top_terms.Pathway,kegg_top_terms.enrichment_score,"output/share/kegg_bar.tex")
### GO enrichment
##add enrichment column
go_top_terms.enrichment_score = broadcast(x->-log(x),go_top_terms.P_DE) 
## remove any underscores from terms
go_top_terms.Term = replace.(go_top_terms.Term,"_"=>" ")

tex_enrichment_bar(go_top_terms.Term,go_top_terms.enrichment_score,"output/share/go_bar.tex")

