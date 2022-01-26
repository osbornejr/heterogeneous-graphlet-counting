using LightGraphs
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


##set up csv string
#csv = "Nodes,Edges,Components,Nodes in largest component,Maximal degree,communities detected\n"
#csv = csv*string(length(vertexlist))*","*string(length(edgelist))*","*string(length(components))*","*string(length(largest...))*","*string(max(degrees...))*","*string(length(unique(community_vertices.group)))*"\n"
#write("$(params.website_dir)/_assets/$(params.page_name)/tableinput/network_stats.csv",csv)

