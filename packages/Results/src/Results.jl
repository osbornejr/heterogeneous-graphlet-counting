module Results
#Preprocessing
using GLMakie,StatsBase
#Network Construction
using Graphs,GraphMakie,NetworkLayout
using ProjectFunctions,GraphletCounting,NetworkConstruction

## Typed representation (TODO check clash between Cairo and GL makie-- make switchable?)
#using MakieTex,CairoMakie,FileIO


function variance_histogram(count_data)
    var_data = vec(var(data_from_dataframe(count_data),dims=2 ))

    #histogram to see where coding and noncoding is represented in variance data
    bins = Int.(floor.(var_data.*10)).+1
    barplot(bins,repeat([1],length(var_data)),stack=bins,color=replace(count_data.transcript_type,"coding"=>1,"noncoding"=>2))
end


function plot_network(adj_matrix;vertex_colors=nothing)
    #define node colors first
    if (vertex_colors == nothing)
        vertex_colors = repeat(["blue"],size(adj_matrix)[1])
    end

    g = Graph(adj_matrix)
    set_theme!(backgroundcolor="#212121",textcolor=:white)
    fig,scene,p = graphplot(g;
                            layout=Spring(dim=3,seed=5),#Spectral(dim=3),
    node_color = vertex_colors,
    node_size = 10,
    edge_color = :white,
    edge_width = .05,
    figure = (size = (1500, 800),)

    )
    scene.show_axis =false
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

#end module
end
