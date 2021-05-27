using Printf
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

function html_table_maker(dataframe::DataFrame,outfile::String)
	#begin table
	#create header from dataframe names
	table = "<table>\n<thead>\n<tr>\n"
	for n in names(dataframe)
		table = table*"<th>$n</th>\n"
	end
	table = table*"</tr>\n</thead>\n"
	##body from rows
	table = table*"<tbody>\n"
	for r in eachrow(dataframe)

	table = table*"<tr>\n"
		for d in r
			if (typeof(d)<:AbstractFloat)
				table = table*"<td>$(@sprintf "%.5g" d)</td>\n"
			else
				table = table*"<td>$(d)</td>\n"
			end
		end
	table = table*"</tr>\n"
	end
	table = table*"</tbody>\n</table>"
	write(outfile,table)


end

