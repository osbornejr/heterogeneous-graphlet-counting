using ProjectFunctions,CSV

##load in information
kegg_top_terms,go_top_terms = ProjectFunctions.biological_validation()

function tex_enrichment_bar(top_names::Vector{String},,enrichment_scores::Vector{Float},outfile::String)

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
    global tex*="($(enrichment_scores[i]),$(top_names[i])) "
end
tex*="};\n\\end{axis}\n\\end{tikzpicture}"
write(outfile,tex)
end

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

go_top_terms.enrichment_score = broadcast(x->-log(x),go_top_terms.P_DE) 

tex_enrichment_bar(go_top_terms.Pathway,go_top_terms.enrichment_score,"output/share/go_bar.tex")
