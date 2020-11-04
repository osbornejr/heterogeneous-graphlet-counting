function cytoscape_elements(vertices::Array{String,2},edges::Array{Pair},output_path::String)
	io = open(output_path, "w")
	println(io, "var Elements = {")
	println(io, "	nodes: [")
	colour = Dict{String,String}("coding"=>"red","noncoding"=>"green")
	for (i,vertex) in enumerate(eachrow(vertices))
		println(io, "		{ data: { id: ",i,",transcript: '",vertex[1],"', type: '",vertex[2],"', colour: '",colour[vertex[2]],"' } },")
	end
	println(io,"	],")
	println(io,"	edges: [")
	for edge in edges 
		println(io, "		{ data: { id: '",string(edge),"', source: ",edge[1],", target: ",edge[2]," } },")
	end
	println(io, "	]")
	println(io, "};")
	close(io)
end

	
