module Results
#Preprocessing
using GLMakie,StatsBase
#Network Construction
using Graphs,GraphMakie,NetworkLayout
using ProjectFunctions,GraphletCounting


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
                            layout=Spectral(dim=3),
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
    ccd = (vertex_typelist.=="coding").* cd
    ncd = (vertex_typelist.=="noncoding").* cd
    cnd = (vertex_typelist.=="coding").* nd
    nnd = (vertex_typelist.=="noncoding").* nd
    f = Figure()

    barplot(f[1,1],dd,repeat([1],length(dd)),stack=dd,color=(vertex_typelist.=="coding").+1)
    hist(f[2,1],cd, bins = 1623)
    hist(f[2,2],nd, bins = 1623)
    hist(f[3,1],ccd, bins = 1623)
    hist(f[3,2],ncd, bins = 1623)
    hist(f[4,1],cnd, bins = 1623)
    hist(f[4,2],nnd, bins = 1623)
    f 

end

#end module
end
