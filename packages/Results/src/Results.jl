module Results
#Preprocessing
using GLMakie,StatsBase
#Network Construction
using Graphs,GraphMakie,NetworkLayout
using ProjectFunctions


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


#end module
end
