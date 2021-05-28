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

function html_table_maker(dataframe::DataFrame,imgs::Array{String,1},outfile::String)
	#begin table container with scoped styling
	table ="""<div id="table container">\n\t<style type="text/css" scoped>\n\t\t.row {\n\t\t\tmargin-left:-5px;\n\t\t\tmargin-right:-5px;\n\t\t}\n\t\t.column {\n\t\t\tfloat: left;\n\t\t\twidth: 50%;\n\t\t}\n\t\ttd {\n\t\t\ttext-align:center;\n\t\t\tvertical-align:middle;\n\t\t}\n\t</style>\n\t<div class="row">\n\t\t<div class="column">\n"""
	for (i,r) in enumerate(eachrow(dataframe))
		##start new column after every 5 rows
		if (mod(i,5)==1)
			#close previous column
			if(i>1)
				table = table*"""\t\t\t\t</tbody>\n\t\t\t</table>\n\t\t</div>\n\t\t<div class="column">\n"""		 
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
				table = table*"\t\t\t\t\t\t<td>$(@sprintf "%.5g" d)</td>\n"
			# add image to table (asssumes image is located at figs/$d.svg)
			elseif (d in imgs)
				table = table*"""\t\t\t\t\t\t<td><img  width=80% src="./figs/$d.svg" title="$d"/></td>\n"""				
				
			else
				table = table*"\t\t\t\t\t\t<td>$(d)</td>\n"
			end
		end
		table = table*"\t\t\t\t\t</tr>\n"
	end
	#close last table
	table = table*"\t\t\t\t</tbody>\n\t\t\t</table>\n"
 	#close columns, rows and container
	table = table*"\t\t</div>\n\t</div>\n</div>"
	write(outfile,table)
end

