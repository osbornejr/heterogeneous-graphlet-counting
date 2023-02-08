using Printf, Luxor, Colors
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


function tex_table_maker(dataframe::DataFrame,outfile::String)
    table = "\\begin{tabular}"


    #number of rows and columns
    rows = size(dataframe)[1]
    columns = size(dataframe)[2]
    
    table *= "{@{}$(repeat("l",columns))@{}}\n"
    table *= "\\toprule\n"

    #add in column names
    for (i,n) in enumerate(names(dataframe))
        #latex doesn't like _
        cn = replace(n,"_"=>"-")
        table *= "$(cn)"
        if (i<columns)
            table *=" & "
        end
    end
    table *= " \\\\\n"
    table *= "\\midrule\n"

    for (i,r) in enumerate(eachrow(dataframe))
        for (j,d) in enumerate(r)
            #format floats and ints
            if (typeof(d)<:AbstractFloat)
                if (d>10000)
                    table *= " $(Int(round(d))) "
                else
                    table *= "$(round(d,digits=4)) "
                end
            else
                table *= "$(d) "
            end 
            if (j<columns)
                table *=" & "
            end
        end
        table *= "\\\\\n"
    end
    table*= "\\bottomrule\n"

    table*= "\\end{tabular}\n"

    write(outfile,table)
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

function threejs_plot(adj_matrix::AbstractArray,vertex_names::Vector{<:AbstractString},colour_by::Vector{<:AbstractString},plot_prefix::String="",colours::Vector{<:AbstractString}=["blue","green"])
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
    plot = graphjs(g,bg = "white")
    saveWidget(plot,paste0(plot_prefix,"_threejs_plot.html"))
    """
end

function tex_boxplot(data::DataFrame,points::Array{Float64,1},out_file::String,out_format::String;ylabel::String="value")
    # generates a modified version of plot discussed at https://tex.stackexchange.com/a/495207
    #initialise tex string
    tex = ""
    if (out_format == "standalone")
        #include standalone preamble
        tex *= "\\documentclass[crop=false]{standalone}\n\\usepackage{pgfplotstable}\n\\usepgfplotslibrary{colorbrewer}\n%\\pgfplotsset{compat=1.16}\n\\usepgfplotslibrary{statistics}\n\n\\begin{document}\n"
    end

    if (out_format in ["standalone","input"])
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
        tex = chop(tex,tail=1)
        tex *= "},\nxtick={"
        for a in 1:length(data[1])
            tex *= string(a)
            if(a<length(data[1]))
               tex*= ","
           end
        end
        tex *= "},\nx tick label style={scale=0.5,font=\\bfseries, rotate=60,,align=center,anchor = east},\nylabel={$ylabel},cycle list/Set3 ]\n\\pgfplotstablegetrowsof{\\datatable}\n\\pgfmathtruncatemacro{\\rownumber}{\\pgfplotsretval-1}\n\\pgfplotsinvokeforeach{0,...,\\rownumber}{ \n\\pgfplotstablegetelem{#1}{min}\\of\\datatable \n\\edef\\mymin{\\pgfplotsretval} \n \n\\pgfplotstablegetelem{#1}{q25}\\of\\datatable \n\\edef\\myql{\\pgfplotsretval} \n \n\\pgfplotstablegetelem{#1}{median}\\of\\datatable \n\\edef\\mymedian{\\pgfplotsretval} \n \n\\pgfplotstablegetelem{#1}{q75}\\of\\datatable \n\\edef\\myqu{\\pgfplotsretval} \n \n\\pgfplotstablegetelem{#1}{max}\\of\\datatable \n\\edef\\mymax{\\pgfplotsretval} \n \n\\typeout{\\mymin,\\myql,\\mymedian,\\myqu,\\mymax} \n\\pgfmathsetmacro{\\mylowerq}{\\myql} \n\\pgfmathsetmacro{\\myupperq}{\\myqu} \n\\edef\\temp{\\noexpand\\addplot+[, \nboxplot prepared={ \n     lower whisker=\\mymin, \n     upper whisker=\\mymax, \n     lower quartile=\\mylowerq, \n     upper quartile=\\myupperq, \n     median=\\mymedian, \n     every box/.style={solid,fill,opacity=0.5}, \n     every whisker/.style={solid }, \n     every median/.style={solid}, \n     }, \n]coordinates {};} \n\\temp \n}\n\\addplot [only marks, mark=o,mark size = 1pt] coordinates{"
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

function tex_merged_boxplot(data_array::Array{DataFrame,1},out_file::String,out_format::String;ylabel::String="value")
    #initialise tex string
    tex = ""
    if (out_format == "standalone")
        #include standalone preamble
        tex *= "\\documentclass[crop=false]{standalone}\n\\usepackage{pgfplotstable}\n\\usepgfplotslibrary{colorbrewer}\n%\\pgfplotsset{compat=1.16}\n\\usepgfplotslibrary{statistics}\n\n\\begin{document}\n"
    end
    ## merge data
    merged_data = vcat(data_array...)
    #get global min and max of all datapoints
    ymin = min(Array(merged_data[!,2:end])...)-0.5
    ymax = max(Array(merged_data[!,2:end])...)+0.5

    if (out_format in ["standalone","input"])
        tex *= "\\begin{tikzpicture}\n\\pgfplotsset{colormap/Set2}\n"
        tex *= "\\begin{axis}[xmin = 0, xmax = $(length(merged_data[1])+1), enlarge x limits ={0},ymin = $ymin, ymax = $ymax, boxplot/draw direction=y,\nxticklabels={"
        #get row names (from column one). IMPORTANTLY, latex cant handle underscores in name. Plots also present better with short labels, so we initialise each typed graphlet
        for a in merged_data[1]
            tex *= replace(replace(replace(chop(split(a,"-")[1],tail=2),"_"=>"-"),"oding"=>""),"on"=>"")*","
        end
        tex = chop(tex,tail=1)
        tex *= "},\nxtick={"
        for a in 1:length(merged_data[1])
            tex *= string(a)*","
        end
        tex = chop(tex,tail=1)
        tex *= "},\nx tick label style={scale=0.5,font=\\bfseries, rotate=60,,align=center,anchor = east},\nylabel={$ylabel},cycle list/Set3 ]\n"
        #cycle through graphlets
        #keep global counts
        tal = 0
        hom_count = 0
        grey_flag = 0
        for data in data_array
            box_data = data[!,Not(:values)]
            #define data for these boxplots
            tex *= "\\pgfplotstableread{%\nx min q25 median q75 max\n"
            for row in eachrow(box_data)
                arr = vec(convert(Array,row))
                for a in arr
                    tex *= string(a)*" "
                end
                tex *= "\n"
            end
            tex *= "}\\datatable\n"
            ##define grey box area every second graphlet
            if ((grey_flag == 1) && (hom_count+1 == length(data_array)))
                tex *= "\\addplot [draw=gray,fill=gray,opacity=0.1] coordinates {($(tal+0.5),$(ymin-1)) ($(tal+0.5),$(ymax+1)) ($(tal+length(data[1])+1),$(ymax+1)) ($(tal+length(data[1])+1),$(ymin-1))};\n"
            elseif (grey_flag == 1)
                tex *= "\\addplot [draw=gray,fill=gray,opacity=0.1] coordinates {($(tal+0.5),$(ymin-1)) ($(tal+0.5),$(ymax+1)) ($(tal+length(data[1])+0.5),$(ymax+1)) ($(tal+length(data[1])+0.5),$(ymin-1))};\n"
                grey_flag = 0
            else
                grey_flag = 1
            end
            ##box plots for this graphlet
            tex*= "\\pgfplotstablegetrowsof{\\datatable}\n\\pgfmathtruncatemacro{\\rownumber}{\\pgfplotsretval-1}\n\\pgfplotsinvokeforeach{0,...,\\rownumber}{ \n\\pgfplotstablegetelem{#1}{min}\\of\\datatable \n\\edef\\mymin{\\pgfplotsretval} \n \n\\pgfplotstablegetelem{#1}{q25}\\of\\datatable \n\\edef\\myql{\\pgfplotsretval} \n \n\\pgfplotstablegetelem{#1}{median}\\of\\datatable \n\\edef\\mymedian{\\pgfplotsretval} \n \n\\pgfplotstablegetelem{#1}{q75}\\of\\datatable \n\\edef\\myqu{\\pgfplotsretval} \n \n\\pgfplotstablegetelem{#1}{max}\\of\\datatable \n\\edef\\mymax{\\pgfplotsretval} \n \n\\typeout{\\mymin,\\myql,\\mymedian,\\myqu,\\mymax} \n\\pgfmathsetmacro{\\mylowerq}{\\myql} \n\\pgfmathsetmacro{\\myupperq}{\\myqu} \n\\edef\\temp{\\noexpand\\addplot+[index of colormap={$(hom_count) of Set2}, \nboxplot prepared={ \n     lower whisker=\\mymin, \n     upper whisker=\\mymax, \n     lower quartile=\\mylowerq, \n     upper quartile=\\myupperq, \n     median=\\mymedian, \n     every box/.style={solid,fill,opacity=0.5}, \n     every whisker/.style={solid }, \n     every median/.style={solid}, \n     }, \n]coordinates {};} \n\\temp \n}\n\\addplot [only marks, mark=o,mark size = 1pt] coordinates{"
            points = data.values
            hom_count +=1 
            for p in points
                tal += 1
                tex*= " ($tal,$p) "
            end
            tex*= "};\n"
        end
        tex *= "\\end{axis}\n\\end{tikzpicture}\n"
    end

    if (out_format == "standalone")
        tex *= "\\end{document}"
    end


    write(out_file,tex)
end 



function graphlet_adjacency_to_edgelist_array(adj::Matrix{Int64})
    return graphlet_adjacency_to_edgelist_array(BitMatrix(adj))
end

function graphlet_adjacency_to_edgelist_array(adj::AbstractMatrix{Bool})
    return [adj[i,j] for i in 1:size(adj)[1]-1 for j in i+1:size(adj)[2]]
end

function graphlet_edgelist_array_to_adjacency(vecc::AbstractVector{Int})
    return graphlet_edgelist_array_to_adjacency(BitVector(vecc))
end

function graphlet_edgelist_array_to_adjacency(vecc::AbstractVector{Bool})
    ## provided vector must be of length equal to triangular number associated with dims of adj/graphlet. Find here and check that length corresponds to a dim that is whole number

    ## inverse triangular number
    m = length(vecc)
    n = sqrt(2*m+0.25)-0.5 
    if !(isinteger(n))
        throw(ArgumentError("provided vector does not correspond to any graphlet dimension n>1. Vector must be of length T(n-1)")) 
    else
        #now we find graphlet dimension
        n = Int(n)+1
    end

    #now add zeros appropriately to vector and then reshape
    adj = []

    for i in 1:(n-1)
        step = n-i
        start = m - Int((step*(step+1))/2) 
        adj = vcat(adj,zeros(Int,i),vecc[start+1:start+step])
    end
    adj = vcat(adj,zeros(Int,n))
    adj = reshape(adj,n,n)
    adj = BitMatrix(adj + adj')
    return adj
end

function generate_heterogeneous_graphlet_list(vecc::AbstractVector{Int},types::Vector{String})
    ##allows for simplifying entry for each edge only. We convert here to a bit matrix and then feed to main function
    #first convert int to Bitvector (all nonzero values converted to true)
    vec = .!(iszero.(vecc))
    adj = graphlet_edgelist_array_to_adjacency(vecc)
    return generate_heterogeneous_graphlet_list(adj,types)
end

function generate_heterogeneous_graphlet_list(vecc::AbstractVector{Bool},types::Vector{String})
    ##allows for simplifying entry for each edge only. We convert here to a bit matrix and then feed to main function
    adj = graphlet_edgelist_array_to_adjacency(vecc)
    return generate_heterogeneous_graphlet_list(adj,types)



end

function generate_heterogeneous_graphlet_list(adj::Matrix{Int64},types::Vector{String})
    return generate_heterogeneous_graphlet_list(BitMatrix(adj),types)
end

function generate_heterogeneous_graphlet_list(adj::BitMatrix,types::Vector{String})
    #method to find all possible permutations of a heterogeneous graphlet given an orbit classification and a set of types.

    #check to make sure types are sorted
    types = sort(types)

    ##deduce orbits from adj (sum across one dim)
    orbits = vec(sum(adj,dims=1)) 

    #generate all possible permutations of the type list
    n = length(types)
    m = length(orbits)

    comb = []
    for i in 1:m
        push!(comb,repeat(vcat([repeat([types[x]],n^(m-i)) for x in 1:n]...),n^(i-1)))
    end
    candidates = hcat(comb...)


    ##find typed orbits for each candidate
    typed_orbits = map(y->map(x->countmap(candidates[y,:][adj[x,:]]),1:m),1:size(candidates,1))    
    ## convert to dict to ignore order and find matching candidates, then finding uniques
    to_dict = countmap.(typed_orbits)
    un_to = unique(to_dict)     

    #for each unique typed orbit, match to first occurence in candidates
    # and add to graphlet list
    graphlets = String[]
    for o in un_to
        c = candidates[findfirst(isequal(o),to_dict),:]
        push!(graphlets,join(c,"_"))
    end
#OUTDATED    #now check each candidate to see if it has a unique symmetry under the orbit structure
#    #initialise empty graphlet list
#    graphlets = String[]
#    for c in eachrow(candidates)
#        flag = 0
#        for o in unique(orbits)
#            test = c[orbits.==o]
#            if !(test == sort(test))
#                #reject as soon as one orbit is not sorted correctly 
#                flag = 1
#
#            elseif(length(unique(test))>1)
#                #if orbit is sorted correctly on two
#            end
#        end
#        if (flag == 0)
#            ##candidate checked out for all orbits, so we add it.
#            push!(graphlets,join(c,"_"))
#        else
#            #reflect by first orbit first, sorting all nodes in that orbit and adjacent ones if necessary   
#
#        end
#
#    end

    return graphlets
end


function draw_tex_graphlet(graphlet_name::String;split_char::String="_",kwargs...)
    slice = string.(split(graphlet_name,split_char)) 
    return draw_tex_graphlet(slice[1:end-1],slice[end];kwargs...)
end

function draw_tex_graphlet(node_schematic::Array{String,1},edge_name::String;kwargs...)
    if (edge_name == "2-path")
        if(length(node_schematic)==0)
            #homogeneous case
            node_schematic = ["one","one"]
        end
        return draw_tex_graphlet(node_schematic,[true];kwargs...)
    elseif (edge_name == "3-path")
        if(length(node_schematic)==0)
            #homogeneous case
            node_schematic = ["one","one","one"]
        end
        return draw_tex_graphlet(node_schematic,[true,false,true];kwargs...)
    elseif (edge_name == "3-tri")
        if(length(node_schematic)==0)
            #homogeneous case
            node_schematic = ["one","one","one"]
        end
        return draw_tex_graphlet(node_schematic,[true,true,true];kwargs...)
    elseif (edge_name == "3-clique")
        if(length(node_schematic)==0)
            #homogeneous case
            node_schematic = ["one","one","one"]
        end
        return draw_tex_graphlet(node_schematic,[true,true,true];kwargs...)
    elseif (edge_name == "4-path")
        if(length(node_schematic)==0)
            #homogeneous case
            node_schematic = ["one","one","one","one"]
        end
        return draw_tex_graphlet(node_schematic,[true,false,true,false,false,true];kwargs...)
    elseif (edge_name == "4-star")
        if(length(node_schematic)==0)
            #homogeneous case
            node_schematic = ["one","one","one","one"]
        end
        return draw_tex_graphlet(node_schematic,[false,true,true,false,false,true];kwargs...)
    elseif (edge_name == "4-tail")
        if(length(node_schematic)==0)
            #homogeneous case
            node_schematic = ["one","one","one","one"]
        end
        return draw_tex_graphlet(node_schematic,[true,true,true,false,false,true];kwargs...)
    elseif (edge_name == "4-cycle")
        if(length(node_schematic)==0)
            #homogeneous case
            node_schematic = ["one","one","one","one"]
        end
        return draw_tex_graphlet(node_schematic,[true,false,true,true,false,true];kwargs...)
    elseif (edge_name == "4-chord")
        if(length(node_schematic)==0)
            #homogeneous case
            node_schematic = ["one","one","one","one"]
        end
        return draw_tex_graphlet(node_schematic,[true,true,true,false,true,true];kwargs...)
    elseif (edge_name == "4-clique")
        if(length(node_schematic)==0)
            #homogeneous case
            node_schematic = ["one","one","one","one"]
        end
        return draw_tex_graphlet(node_schematic,[true,true,true,true,true,true];kwargs...)
    else
        throw(ArgumentError("$edge_name not recognised as a valid default schematic. Please provide edge schematic explicitly"))
    end
end

function draw_tex_graphlet(node_schematic::Vector{String},edge_schematic::AbstractVector{Bool};out_file::String,colours::Vector{String}=["black"])
    #function to create tikz drawings of a given graphlet
   
    #initialise tex text
    tex = "\\begin{tikzpicture}[main_node/.style={circle,fill=black,minimum size=2em,inner sep=3pt]}]\n"

    #NODES
    order = length(node_schematic)
    node_names = sort(unique(node_schematic))
    ##colour mapping
    #first check that there are enough provided colours. If not, default to black
    if(length(node_names)>length(colours))     
        @info "Insufficient number of colours provided, defaulting to black for all nodes"
        colours = fill("black",length(node_names))
    end
    colour_schematic = copy(node_schematic)
    for (i,node) in enumerate(node_names)
        replace!(colour_schematic,node=>colours[i])  
    end

    ##add in node points (with colours)
    tex = tex*"\\node[main_node,fill=$(colour_schematic[1])] (1) at (0,0) {};\n"

    tex = tex*"\\node[main_node,fill=$(colour_schematic[2])] (2) at (0, -1.5)  {};\n"
    
    if (order>2)
        tex = tex*"\\node[main_node,fill=$(colour_schematic[3])] (3) at (1.5, -1.5) {};\n"
    end
    if (order>3)
        tex = tex*"\\node[main_node,fill=$(colour_schematic[4])] (4) at (1.5, 0) {};\n"
    end
    #add dummy node for 2-path to keep size fitting
    if (order==2)
        tex = tex*"\\node[main_node,fill=white,opacity=0] (3) at (1.5, -1.5) {};\n"
    end
    


    #EDGES
    #check edge schematic is valid for order
    s_l = Int(order*(order-1)/2)
    if (length(edge_schematic)!= s_l)
        throw(DimensionMismatch("Graphlet edge schematic is not valid. For order $order graphlet, please provide edge schematic of length $s_l."))
    end
    if (edge_schematic[1])
        tex = tex*"\\draw (1) -- (2);\n"
    end
    if (order>2)
        if (edge_schematic[2])
            tex = tex*"\\draw (1) -- (3);\n"
        end
        if (edge_schematic[3])
            tex = tex*"\\draw (2) -- (3);\n"
        end
    end
    if (order>3)
        if (edge_schematic[4])
            tex = tex*"\\draw (1) -- (4);\n"
        end
        if (edge_schematic[5])
            tex = tex*"\\draw (2) -- (4);\n"
        end
        if (edge_schematic[6])
            tex = tex*"\\draw (3) -- (4);\n"
        end
    end
        
    tex = tex*"\\end{tikzpicture}\n"    
    write(out_file,tex)
end


function draw_graphlet(adj::AbstractMatrix;kwargs...)
    return draw_graphlet(["v" for i in 1:size(adj)[1]], graphlet_adjacency_to_edgelist_array(adj);kwargs...)
end

function draw_graphlet(node_schematic::Vector{String},adj::AbstractMatrix;kwargs...)
    return draw_graphlet(node_schematic,graphlet_adjacency_to_edgelist_array(adj);kwargs...)
end

function draw_graphlet(graphlet_name::String;split_char::String="_",kwargs...)
    slice = string.(split(graphlet_name,split_char)) 
    
    if (length(slice)==1)
        edge_name = slice[1]
        if (edge_name == "2-path")
            return draw_graphlet(["v","v"],[true];kwargs...)
        elseif (edge_name == "3-path")
            return draw_graphlet(["v","v","v"],[true,false,true];kwargs...)
        elseif (edge_name == "3-tri")
            return draw_graphlet(["v","v","v"],[true,true,true];kwargs...)
        elseif (edge_name == "3-clique")
            return draw_graphlet(["v","v","v"],[true,true,true];kwargs...)
        elseif (edge_name == "4-path")
            return draw_graphlet(["v","v","v","v"],[true,false,false,true,false,true];kwargs...)
        elseif (edge_name == "4-star")
            return draw_graphlet(["v","v","v","v"],[false,true,false,true,false,true];kwargs...)
        elseif (edge_name == "4-tail")
            return draw_graphlet(["v","v","v","v"],[true,true,false,true,false,true];kwargs...)
        elseif (edge_name == "4-cycle")
            return draw_graphlet(["v","v","v","v"],[true,false,true,true,false,true];kwargs...)
        elseif (edge_name == "4-chord")
            return draw_graphlet(["v","v","v","v"],[true,true,true,true,false,true];kwargs...)
        elseif (edge_name == "4-clique")
            return draw_graphlet(["v","v","v","v"],[true,true,true,true,true,true];kwargs...)
        else
            throw(ArgumentError("$edge_name not recognised as a valid default schematic. Please provide edge schematic explicitly"))
        end
    else
        return draw_graphlet(slice[1:end-1],slice[end];kwargs...)
    end
end

function draw_graphlet(node_schematic::Vector{String},edge_name::String;kwargs...)
    if (edge_name == "2-path")
        return draw_graphlet(node_schematic,[true];kwargs...)
    elseif (edge_name == "3-path")
        return draw_graphlet(node_schematic,[true,false,true];kwargs...)
    elseif (edge_name == "3-tri")
        return draw_graphlet(node_schematic,[true,true,true];kwargs...)
    elseif (edge_name == "3-clique")
        return draw_graphlet(node_schematic,[true,true,true];kwargs...)
    elseif (edge_name == "4-path")
        return draw_graphlet(node_schematic,[true,false,false,true,false,true];kwargs...)
    elseif (edge_name == "4-star")
        return draw_graphlet(node_schematic,[true,false,false,true,true,false];kwargs...)
    elseif (edge_name == "4-tail")
        return draw_graphlet(node_schematic,[true,false,false,true,true,true];kwargs...)
    elseif (edge_name == "4-cycle")
        return draw_graphlet(node_schematic,[true,false,true,true,false,true];kwargs...)
    elseif (edge_name == "4-chord")
        return draw_graphlet(node_schematic,[true,true,false,true,true,true];kwargs...)
    elseif (edge_name == "4-clique")
        return draw_graphlet(node_schematic,[true,true,true,true,true,true];kwargs...)
    else
        throw(ArgumentError("$edge_name not recognised as a valid default schematic. Please provide edge schematic explicitly"))
    end
end

function draw_graphlet(node_schematic::Vector{String},edge_schematic::Vector{Int};kwargs...)
    return draw_graphlet(node_schematic,BitVector(edge_schematic);kwargs...)
end

function draw_graphlet(node_schematic::Vector{String},edge_schematic::AbstractVector{Bool};dim::Int=50,rotation::Float64=0.0,node_colours::Array{String,1}=["hotpink","gold","skyblue","limegreen"],line_colour::String = "lightgrey",file::Union{Symbol,String}=:svg) 
 #function to create graphlet images programatically using Luxor tools.
 #`file` can be either a filepath string to save a hard copy of image, or a symbol (either :svg or :png) to only create image in memory (useful for cases (i.e. Pluto) where separate file artefacts are a drawback.  
    Drawing(dim,dim,file) 
    origin()
    order = length(node_schematic)
    #rotation: can be set manually, or otherwise use defaults for each order (up to order 4 atm).
    if (rotation!=0.0)
        r = rotation
    elseif (order == 2)
        r = -pi/4
    elseif (order == 3)
        r = -3*pi/2
    elseif (order == 4)
        r = 3*pi/4
    else
        r = 0.0
    end

    

    Luxor.rotate(r)
   
    #create node points around origin
    corners = Luxor.ngon(Point(0, 0), dim/3,order, vertices=true)
    polyreflect!(corners)    

    ##EDGES
    #draw lines according to edge schematic 
    setcolor(line_colour)
    #check edge schematic is valid for order
    s_l = Int(order*(order-1)/2)
    if (length(edge_schematic)!= s_l)
        throw(DimensionMismatch("Graphlet edge schematic is not valid. For order $order graphlet, please provide edge schematic of length $s_l."))
    end
    #count edges for edge schematic
    s_c = 0
    for i in 1:length(corners)
        for j in (i+1):length(corners)
            s_c += 1
            if (edge_schematic[s_c])
                Luxor.line(corners[i],corners[j], action=:stroke)
            end
        end
    end
    #NODES
    #get sorted list of types
    types = sort(unique(node_schematic))
    # set palette to appropriate length, ideally based on node_colours
    if (length(types)<=length(node_colours))
        col_pal = node_colours[1:length(types)]
    else
        #if more colors are required than provided, create new palette to ensure distinguishable colours
        col_pal = distinguishable_colors(length(types)) 
    end

    #now add chosen colours to node schematic in place of types
    #first make sure colours are in same order as types (to avoid erroneous swapping)
    sort!(col_pal)
    for (i,t) in enumerate(types)
        replace!(node_schematic,t=>col_pal[i])
    end
    
    #check that length of node_schematic matches order
    if (length(node_schematic)!=order)
        @error "Graphlet node schematic is not valid. Ensure length of node schematic is equal to $order."
    end
    
    
    #draw and colour nodes according to node schematic
    for i in 1:length(corners)
        setcolor(node_schematic[i])
        Luxor.circle.(corners[i], dim/10, action=:fill)
    end
    finish()
    preview()
end

function three_js_network(adj_matrix::AbstractArray,vertex_names::Vector{<:AbstractString},colours::Vector{<:AbstractString})

    R"""
    sapply(names(sessionInfo()$otherPkgs),function(pkg) detach(paste0('package:',pkg),character.only =T,force = T));rm(list=ls())
    """
    @rput adj_matrix
    @rput vertex_names
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
    
    g <- graph_from_data_frame(edges,directed = F, vertices)
    #
    ##delete zero degree vertices
    g <- delete.vertices(g,igraph::degree(g)==0)
    ##find maximum connected component
    g <- decompose(g, mode = "strong", max.comps = NA, min.vertices = 10)[[1]]
    #
    ##add colours to graph
    vertices <- vertices %>% mutate(group = as_factor(colours),color = group)
    colour_palette = brewer.pal(name = "Spectral", n = length(unique(colours)))
    levels(vertices$color) <- colour_palette
    edges <- as_tibble(as.data.frame(get.edgelist(g)))
    g <- graph_from_data_frame(edges,directed = F, vertices)
    vertex_attr(g,"size") <- 0.5
    plot = graphjs(g,bg = "white")
    """
end
