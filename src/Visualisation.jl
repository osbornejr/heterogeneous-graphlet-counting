using Printf, LaTeXStrings
function cytoscape_elements(vertices::Array{String,2},edges::Array{Pair},output_path::String)
    io = open(output_path, "w")
    println(io, "var Elements = {")
    println(io, "   nodes: [")
    colour = Dict{String,String}("coding"=>"red","noncoding"=>"green")
    for (i,vertex) in enumerate(eachrow(vertices))
        println(io, "       { data: { id: ",i,",transcript: '",vertex[1],"', type: '",vertex[2],"', colour: '",colour[vertex[2]],"' } },")
    end
    println(io,"    ],")
    println(io,"    edges: [")
    for edge in edges 
        println(io, "       { data: { id: '",string(edge),"', source: ",edge[1],", target: ",edge[2]," } },")
    end
    println(io, "   ]")
    println(io, "};")
    close(io)
end

function html_table_maker(dataframe::DataFrame,outfile::String;imgs::Array{String,1}=[],figpath::String=".figs/")
    #begin table container with scoped styling
    table ="""<div id="table container">\n\t<style type="text/css" scoped>\n\t\t.row {\n\t\t\tmargin-left:-5px;\n\t\t\tmargin-right:-5px;\n\t\t}\n\t\t.column {\n\t\t\tfloat: left;\n\t\t\twidth: 50%;\n\t\t}\n\t\ttd {\n\t\t\ttext-align:center;\n\t\t\tvertical-align:middle;\n\t\t}\n\t</style>\n\t<div class="row">\n\t\t<div class="column">\n"""
    for (i,r) in enumerate(eachrow(dataframe))
        ##start new column after every 5 rows
        if (mod(i,5)==1)
            #close previous column
            if(i>1)
                #close table
                table = table*"""\t\t\t\t</tbody>\n\t\t\t</table>\n"""       
                #add legend to bottom of first row
                table = table*"""\t\t<div class="column">\n\t\t\t<table align="left">\n\t\t\t\t<tr>\n\t\t\t\t\t<td bgcolor="#fc3c7a">&nbsp;&nbsp;&nbsp;</td>\n\t\t\t\t\t<td>coding</td>\n\t\t\t\t\t<td bgcolor="#f4cd16">&nbsp;&nbsp;&nbsp;</td>\n\t\t\t\t\t<td>noncoding</td>\n\t\t\t\t</tr>\n\t\t\t</table>\n\t\t</div>\n"""
                #close column
                table = table*"""\n\t\t</div>\n\t\t<div class="column">\n"""         
            end
            #begin table...create header from dataframe names
            table = table*"\t\t\t<table>\n\t\t\t\t<thead>\n\t\t\t\t\t<tr>\n"
            for n in names(dataframe)
                table = table*"\t\t\t\t\t\t<th>$n</th>\n"
            end
            table = table*"\t\t\t\t\t</tr>\n\t\t\t\t</thead>\n"
            ##body for next 5 rows
            table = table*"\t\t\t\t<tbody>\n"
        end
        ##begin row entry
        table = table*"\t\t\t\t\t<tr>\n"
        for d in r
            # format floats and ints    
            if (typeof(d)<:AbstractFloat)
                table = table*"\t\t\t\t\t\t<td>$(@sprintf "%.3g" d)</td>\n"
            # add image to table (asssumes image is located at figs/$d.svg)
            elseif (d in imgs)
                table = table*"""\t\t\t\t\t\t<td><img  width=80% src="$figpath/$d.svg" title="$d"/></td>\n"""               
                
            else
                table = table*"\t\t\t\t\t\t<td>$(d)</td>\n"
            end
        end
        table = table*"\t\t\t\t\t</tr>\n"
    end
    #close last table
    table = table*"\t\t\t\t</tbody>\n\t\t\t</table>\n"
    #close column 
    table = table*"\t\t</div>\n"
    #close row
    table = table*"\t</div>\n"
    #close container
    table = table*"</div>\n"
    write(outfile,table)
end

function threejs_plot(adj_matrix::AbstractArray,vertex_names::Array{String,1},colour_by::Array{String,1},plot_prefix::String,colours::Array{String,1})
    R"""
    sapply(names(sessionInfo()$otherPkgs),function(pkg) detach(paste0('package:',pkg),character.only =T,force = T));rm(list=ls())
    """
    @rput adj_matrix
    @rput vertex_names
    @rput colour_by
    @rput plot_prefix   
    @rput colours
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
    vertices <- vertices %>% mutate(color = as_factor(colour_by))
    g <- graph_from_data_frame(edges,directed = F, vertices)
    #
    ##delete zero degree vertices
    g <- delete.vertices(g,igraph::degree(g)==0)
    ##find maximum connected component
    #g <- decompose(g, mode = "strong", max.comps = NA, min.vertices = 10)[[1]]
    #
    ##get only connected vertices
    #vertices <- vertices %>% filter(name %in% names(V(g)))
    
    ##add colours to graph
    #colour_palette = brewer.pal(name = "Spectral", n = length(levels(vertices$color)))
    colour_palette = colours
    levels(vertices$color) <- colour_palette
    edges <- as_tibble(as.data.frame(get.edgelist(g)))
    g <- graph_from_data_frame(edges,directed = F, vertices)
    vertex_attr(g,"size") <- 0.5
    plot = graphjs(g,bg = "white");
    saveWidget(plot,paste0(plot_prefix,"_threejs_plot.html"))
    """

    @info "threejs plot saved as HTML page at $(plot_prefix)."
end


function tex_boxplot(data::DataFrame,points::Array{Float64,1},out_file::String,out_format::String,ylabel::String="value")
    if (out_format == "standalone")
        #include standalone preamble
        tex = "\\documentclass[crop=false]{standalone}\n\\usepackage{pgfplotstable}\n\\usepgfplotslibrary{colorbrewer}\n%\\pgfplotsset{compat=1.16}\n\\usepgfplotslibrary{statistics}\n\n\\begin{document}\n"
    end

    if (out_format == "standalone" || "input")
        tex *= "\\begin{tikzpicture}\n\\pgfplotstableread{%\nx min q25 median q75 max\n"
         for row in eachrow(data)
             arr = vec(convert(Array,row))
             for a in arr
                 tex *= string(a)*" "
             end
             tex *= "\n"
         end
        tex *= "}\\datatable\n\\begin{axis}[boxplot/draw direction=y,\nxticklabels={"
        #get row names (from column one). IMPORTANTLY, latex cant handle underscores in name. Plots also present better with short labels, so we initialise each typed graphlet
        for a in data[1]
            tex *= replace(replace(replace(a,"_"=>"-"),"oding"=>""),"on"=>"")*","
        end
        tex *= "},\nxtick={"
        for a in 1:length(data[1])
            tex *= string(a)
            if(a<length(data[1]))
               tex*= ","
           end
        end
        tex *= "},\nx tick label style={scale=0.5,font=\\bfseries, rotate=60,,align=center},\nylabel={$ylabel},cycle list/Set3 ]\n\\pgfplotstablegetrowsof{\\datatable}\n\\pgfmathtruncatemacro{\\rownumber}{\\pgfplotsretval-1}\n\\pgfplotsinvokeforeach{0,...,\\rownumber}{ \n\\pgfplotstablegetelem{#1}{min}\\of\\datatable \n\\edef\\mymin{\\pgfplotsretval} \n \n\\pgfplotstablegetelem{#1}{q25}\\of\\datatable \n\\edef\\myql{\\pgfplotsretval} \n \n\\pgfplotstablegetelem{#1}{median}\\of\\datatable \n\\edef\\mymedian{\\pgfplotsretval} \n \n\\pgfplotstablegetelem{#1}{q75}\\of\\datatable \n\\edef\\myqu{\\pgfplotsretval} \n \n\\pgfplotstablegetelem{#1}{max}\\of\\datatable \n\\edef\\mymax{\\pgfplotsretval} \n \n\\typeout{\\mymin,\\myql,\\mymedian,\\myqu,\\mymax} \n\\pgfmathsetmacro{\\mylowerq}{\\myql} \n\\pgfmathsetmacro{\\myupperq}{\\myqu} \n\\edef\\temp{\\noexpand\\addplot+[, \nboxplot prepared={ \n     lower whisker=\\mymin, \n     upper whisker=\\mymax, \n     lower quartile=\\mylowerq, \n     upper quartile=\\myupperq, \n     median=\\mymedian, \n     every box/.style={solid,fill,opacity=0.5}, \n     every whisker/.style={solid }, \n     every median/.style={solid}, \n     }, \n]coordinates {};} \n\\temp \n}\n\\addplot [only marks, mark=o,mark size = 1pt] coordinates{"
        for (i,p) in enumerate(points)
            tex*= " ($i,$p) "
        end
    tex*= "};\n\\end{axis}\n\\end{tikzpicture}\n"

    end


    if (out_format == "standalone")
        tex *= "\\end{document}"
    end


    write(out_file,tex)
end

