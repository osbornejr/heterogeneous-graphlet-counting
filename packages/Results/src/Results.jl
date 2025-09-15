module Results
#Preprocessing
using StatsBase,DataFrames
#Network Construction
using Graphs,GraphMakie,NetworkLayout
using ProjectFunctions,GraphletCounting,NetworkConstruction

## Typed representation (TODO check clash between Cairo and GL makie-- make switchable?)
using MakieTeX,FileIO,Luxor
using CairoMakie
#using GLMakie


function variance_histogram(count_data)
    var_data = vec(var(data_from_dataframe(count_data),dims=2 ))

    #histogram to see where coding and noncoding is represented in variance data
    bins = Int.(floor.(var_data.*10)).+1
    barplot(bins,repeat([1],length(var_data)),stack=bins,color=replace(count_data.transcript_type,"coding"=>:pink,"noncoding"=>:yellow))
end


function plot_network(adj_matrix;layout="Spring",dim=2,engine="Cairo",vertex_colors=nothing,groupby=nothing)
    #define node colors first
    if (vertex_colors == nothing)
        vertex_colors = repeat(["blue"],size(adj_matrix)[1])
    end
    #dim of graph
    d = dim
    
    #build graph object before setting layout method
    g = Graph(adj_matrix)


    #set layout here
    if layout == "Spring"
        lo = Spring(dim=d,seed=10)
    elseif layout == "Spectral"
        lo = Spectral(dim=d)
    elseif layout == "Stress"
        lo = Stress(dim=d)
    elseif layout == "Grouped"
        if groupby == nothing
            throw(ArgumentError("for Grouped layout, a groupby vector must be supplied"))
        else
            lo = community_force_layout(g, groupby)
        end
    else
        throw(ArgumentError("supplied layout not recognise:w
                            d"))
    end


    set_theme!(backgroundcolor="#212121",textcolor=:white)
    fig,scene,p = graphplot(g;
                            #layout=Spring(dim=d,seed=10),#Spectral(dim=d),
                            layout =lo,
                            node_color = vertex_colors,
                            node_size = 10,
                            edge_color = :grey,#:white,
                            edge_width = .01,
                            figure = (size = (1500, 800),)
                           )
    if d == 3
        scene.show_axis =false
    elseif d == 2
        hidedecorations!(scene)
    end

    return fig
end

function add_to_fig(fig)
    ax2 = Axis(fig[1, 2])  # Add another axis to the right of the graphplot
    lines!(ax2, 1:10, rand(10))

    display(fig)
end

function typed_degree_distribution(vertex_typelist,edgelist)
    tdd = GraphletCounting.typed_degree_distribution(vertex_typelist,edgelist)
    ##get raw degree distribution
    dd = sum.(values.(tdd))
    #coding degrees TODO generalise to any types and any number of types
    cd = first.(collect.(values.(tdd)))
    nd = last.(collect.(values.(tdd)))
    #break down further to differentiate type of source
    ccd = cd[(vertex_typelist.=="coding")]
    ncd = cd[(vertex_typelist.=="noncoding")]
    cnd = nd[(vertex_typelist.=="coding")]
    nnd = nd[(vertex_typelist.=="noncoding")]
    f = Figure()
    ax1 = Axis(f[1, 1], title = "Degree Distribution")
    barplot!(ax1,dd,repeat([1],length(dd)),stack=dd,color=(vertex_typelist.=="coding").+1)
    ax2 = Axis(f[2, 1], title = "Coding Degree Distribution")
    barplot!(ax2,cd,repeat([1],length(cd)),stack=cd,color=(vertex_typelist.=="coding").+1)
    ax3 = Axis(f[2, 2], title = "Noncoding Degree Distribution")
    barplot!(ax3,nd,repeat([1],length(nd)),stack=nd,color=(vertex_typelist.=="coding").+1)
    ax4 = Axis(f[3, 1], title = "Coding to coding Degree Distribution")
    barplot!(ax4,ccd,repeat([1],length(ccd)),stack=ccd,color="gold")
    ax5 = Axis(f[3, 2], title = "Coding to noncoding Degree Distribution")
    barplot!(ax5,ncd,repeat([1],length(ncd)),stack=ncd,color="purple")
    ax6 = Axis(f[4, 1], title = "Noncoding to coding Degree Distribution")
    barplot!(ax6,cnd,repeat([1],length(cnd)),stack=cnd,color="gold")
    ax7 = Axis(f[4, 2], title = "Noncoding to noncoding Degree Distribution")
    barplot!(ax7,nnd,repeat([1],length(nnd)),stack=nnd,color="purple")
    f 

end

function typed_representation_results(t_r_output::Vector{DataFrame};colour_mapping=nothing)
    ## if no colour mapping is provided, need to set one here based on all types in output 
    if isnothing(colour_mapping)
        types = sort(String.(unique(vcat(map(x->x[1:end-1],split.(vcat(map(x->x.Graphlet,t_r_output)...),"_"))...))))
        node_colours = ["hotpink","gold","skyblue","limegreen"]
        colour_mapping = Dict(Pair.(types,node_colours[1:length(types)]))
    end

    fig = Figure(
                 backgroundcolor = RGBf(0.90, 0.90, 0.90),
                 size = (700, 700))
    ## re-sort output so that 3 node graphlets are first, and that within each graphlet order they are sorted from least to most typed orbits
    three_node = collect(1:length(t_r_output))[occursin.("3-",vcat(map(x->String.(unique(map(x->x[end],split.(x.Graphlet,"_")))),t_r_output)...))]
    four_node = collect(1:length(t_r_output))[occursin.("4-",vcat(map(x->String.(unique(map(x->x[end],split.(x.Graphlet,"_")))),t_r_output)...))]
    t_r_output = [sort(t_r_output[three_node],by=nrow)...,sort(t_r_output[four_node],by=nrow)...]    
    for (i,df) in enumerate(t_r_output)
        ax = Axis(fig[i,1],
                  xgridvisible=false,
                  ygridvisible=false,
                  xticklabelsvisible = false,
                  yticklabelsvisible = false,
                  xticksvisible = false,
                  yticksvisible = false
                 )  
        hidespines!(ax)
        #sort by graphlet to maintain consistency on each axis
        sort!(df,:Graphlet)
        for (j,g) in enumerate(df.Graphlet)
            #need to set colours here based on graphlet order so that coding and noncoding colouring is consistent.
            
            ##identify overrepresented graphlets
            @drawsvg begin setline(6);sethue("#42f587"); box(O, 150, 150, 10, action = :stroke) end
            if df[j,:p_value]<0.05
                svg= SVGDocument(svgstring())
                scatter!(ax,j,1,marker=Cached(svg), markersize = 200)
            end
            ##identify underrepresented graphlets
            @drawsvg begin setline(6);sethue("#f74859"); box(O, 150, 150, 10, action = :stroke) end
            if df[j,:p_value]>0.95
                svg= SVGDocument(svgstring())
                scatter!(ax,j,1,marker=Cached(svg), markersize = 200)
            end

            NetworkConstruction.draw_graphlet(g,colour_mapping=colour_mapping)
            svg= SVGDocument(svgstring())
            scatter!(ax,j,1,marker=Cached(svg), markersize = 50)
        end
    end
    #scatter(fig[1, 1], 1, 1, marker=Cached(svg), markersize = 50)
    fig

end

function community_force_layout(g::Graph, communities::Vector{Any})
    num_nodes = nv(g)
    positions = Vector{Makie.Point2f0}(undef, num_nodes)
    
    # turn communities into integers if required
    if !(unique(typeof.(unique(communities)))... == Int)
        communities = Int.(replace(communities,[Pair(c,i) for (i,c) in enumerate(unique(communities))]...))
    end

    # Get unique communities and arrange them in a grid
    unique_comms = sort(unique(communities))
    ncomms = length(unique_comms)
    grid_cols = ceil(Int, sqrt(ncomms))
    grid_rows = ceil(Int, ncomms / grid_cols)
    comm_positions = Dict{Int, Makie.Point2f0}()

    spacing = 0.02  # distance between community centers
    for (i, comm) in enumerate(unique_comms)
        row = (i - 1) ÷ grid_cols
        col = (i - 1) % grid_cols
        comm_positions[comm] = Makie.Point2f0(col * spacing, -row * spacing)
    end

    for comm in unique_comms
        # Get nodes in the community
        nodes = findall(x -> x == comm, communities)
        subg = induced_subgraph(g, nodes)
        alg = Spectral(dim=2)
        subpos = alg(subg[1])

        # Re-map subgraph node indices to original graph node indices
        for (i, node) in enumerate(nodes)
            # Offset local layout by community center
            offset = comm_positions[comm]

            positions[node] = subpos[i] .+ offset
        end
    end

    return positions
end



#end module
end
