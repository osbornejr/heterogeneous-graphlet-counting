using ProjectFunctions,CSV

##bio validation table of top kegg terms associated with network
top_terms = ProjectFunctions.biological_validation()[3]
## calculate enrichment score (invert p-values on log scale)
top_terms.enrichment_score = broadcast(x->-log(x),top_terms.P_DE) 
## remove homo sapiens tail from every pathway (mention in caption)
top_terms.Pathway = replace.(top_terms.Pathway," - Homo sapiens (human)"=>"")
#rekmove other longwinded subcats
top_terms.Pathway = replace.(top_terms.Pathway," - multiple diseases"=>"")
top_terms.Pathway = replace.(top_terms.Pathway," - reactive oxygen species"=>"")

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
tex*= join(sort(top_terms,:enrichment_score).Pathway,",")*"},\n"
tex*="ytick=data]\n"
##add coordinates for bars
tex*= "\\addplot coordinates {"
for i in eachrow(sort(top_terms,:enrichment_score))
    global tex*="($(i.enrichment_score),$(i.Pathway)) "
end
tex*="};\n\\end{axis}\n\\end{tikzpicture}"
write("output/share/bio_bar.tex",tex)
    



