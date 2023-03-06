### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 63f79e1a-b0d8-11ed-1669-a770b3844775
# ╠═╡ show_logs = false
begin
cwd = ENV["PWD"];	
	import Pkg
	Pkg.activate(mktempdir())
    Pkg.add([
        Pkg.PackageSpec(name="DataFrames", version="1.3"),
		Pkg.PackageSpec(name="Revise"),
		Pkg.PackageSpec(name="CommonMark"),
		Pkg.PackageSpec(name="PlutoUI"),
		Pkg.PackageSpec(name="DataStructures"),
		Pkg.PackageSpec(name="StatsBase"),
		Pkg.PackageSpec(name="Gadfly"),
		Pkg.PackageSpec(name="LinearAlgebra")
    ])
	using Revise	
	using CommonMark
	using PlutoUI
	using LinearAlgebra
	using DataFrames
	using DataStructures
	using StatsBase
	using Gadfly
	## load dev packages last
	dev = Pkg.develop
	dev(path=cwd*"/packages/DataPreprocessing")
	dev(path=cwd*"/packages/NetworkConstruction")
	dev(path=cwd*"/packages/GraphletCounting")
	dev(path=cwd*"/packages/GraphletAnalysis")
	dev(path=cwd*"/packages/ProjectFunctions")
	using DataPreprocessing
	using NetworkConstruction
	using GraphletCounting
	using GraphletAnalysis
	using ProjectFunctions
end

# ╔═╡ 5de6c623-1c96-41ba-9bda-534bfbeb8928
# ╠═╡ show_logs = false
load_config(cwd*"/config/run-files/GSE68559.yaml")

# ╔═╡ 170965d2-0728-4baf-950e-f2b56739427d
# ╠═╡ show_logs = false
raw_counts,processed_counts,similarity_matrix,adj_matrix,network_counts,vertexlist,edgelist = get_output_data();

# ╔═╡ 02211a8a-87a4-4bd5-acf8-74cd27b902c6
dd = GraphletCounting.typed_degree_distribution(vertexlist,edgelist);

# ╔═╡ 5e6beb8d-47fc-4a1b-926e-ae17a51b2ee4
plot_dd_data = DataFrame(coding = map(x->DefaultDict(0,x)["coding"],dd),noncoding = map(x->DefaultDict(0,x)["noncoding"],dd));

# ╔═╡ db25a7ae-26b0-40da-97d3-06edb70c70d9
bc = 50

# ╔═╡ 53770902-0cee-4880-b480-e475e158a369
# ╠═╡ show_logs = false
plot(
	layer(plot_dd_data[vertexlist.=="coding",:],x = "coding",Geom.line,Stat.histogram(bincount=bc),color=["coding->coding"]),
	layer(plot_dd_data[vertexlist.=="coding",:],x = "noncoding",Geom.line,Stat.histogram(bincount=bc),color=["coding->noncoding"]),
	layer(plot_dd_data[vertexlist.=="noncoding",:],x = "coding",Geom.line,Stat.histogram(bincount=bc),color=["noncoding->coding"]),
	layer(plot_dd_data[vertexlist.=="noncoding",:],x = "noncoding",Geom.line,Stat.histogram(bincount=bc),color=["noncoding->noncoding"]),
	Guide.title("Typed degree distribution"),
	Guide.xlabel("degree"),
	Guide.ylabel("frequency"),
	Theme(grid_line_width = 0mm))

# ╔═╡ 0ed6a565-9abc-428b-bf0f-81be4a0b1f9c
components = NetworkConstruction.network_components(adj_matrix) 

# ╔═╡ 392062f3-e17f-4532-9337-692c229944a1
 processed_data = DataPreprocessing.data_from_dataframe(processed_counts,"data")

# ╔═╡ 96e7e7aa-ec6a-44fd-8914-70c3f44c771f
dodgy_component_data = processed_data[components[2],:]

# ╔═╡ df684d90-6453-41f8-b0b2-7ca47a4f9ce3
#create binned counts to compare in plot (only for component 2)
begin
	binned_counts = processed_counts[components[2],:];
	binned_counts[:,2:99] = NetworkConstruction.discretise(processed_data[components[2],:],nbins = 50);
	binned_long_data = DataFrame(
	names = repeat(binned_counts.transcript_id,length(2:99)),
	var = stack(binned_counts[:,2:99])[:,1],
	val = stack(binned_counts[:,2:99])[:,2]
);
	##plot showing bins across samples for each node in network. Coloured by two components in network
Gadfly.plot(binned_long_data,x=:var,y=:val,group=:names,Geom.line,Guide.xticks(label=false),
	Guide.xlabel("sample"),
	Guide.ylabel("bin"),
	Guide.title("All bin values across samples for nodes in network"),
	Theme(grid_line_width=0mm));
end

# ╔═╡ 05cec146-5a7b-4bc6-8797-d7f89f5b8b8b
comp_col = string.(in.(1:2582,Ref(components[2])).+1);

# ╔═╡ 02d67238-533e-4c36-85c4-f828f130e340
#comp_col[top_count] = "782";

# ╔═╡ 7882c063-a58f-4747-a3b0-54c4a756bcf4
top_count = findall(processed_data.==max(processed_data...))[1][1]

# ╔═╡ 7e828b37-d757-4c23-bf21-aa821e776291
processed_long_data = DataFrame(
	names = repeat(processed_counts.transcript_id,length(2:99)),
	component = repeat(comp_col,length(2:99)),
	var = stack(processed_counts[:,2:99])[:,1],
	val = stack(processed_counts[:,2:99])[:,2]
);

# ╔═╡ f179209c-cb71-437c-95c6-ba39513ee394
##plot showing counts across samples for each node in network. Coloured by two components in network
Gadfly.plot(processed_long_data,x=:var,y=:val,group=:names,color= :component,Geom.line,Guide.xticks(label=false),
	Guide.xlabel("sample"),
	Guide.ylabel("count"),
	Guide.title("All processed counts in network"),
	Theme(grid_line_width=0mm))


# ╔═╡ 69a0ad82-e6e7-4832-b393-95ff5561722e
processed_long_data_comp_2 = DataFrame(
	names = repeat(processed_counts.transcript_id[components[2]],length(2:99)),
	var = stack(processed_counts[components[2],2:99])[:,1],
	val = stack(processed_counts[components[2],2:99])[:,2]
);

# ╔═╡ c97ea9f4-c179-4b8e-9abe-d180b75e1abb
##separate plot for 2nd component (if needed)
Gadfly.plot(processed_long_data_comp_2,x=:var,y=:val,group=:names,Geom.line,Guide.xticks(label=false),
	Guide.xlabel("sample"),
	Guide.ylabel("count"),
	Guide.title("2nd component (clique) processed counts"),
	Theme(key_position = :none,grid_line_width=0mm));

# ╔═╡ d651be50-9d95-49fc-bc70-dfc119947a6c
map_component_back_to_raw = findall(x->x in processed_counts.transcript_id[components[2]],raw_counts.transcript_id);

# ╔═╡ a854f969-4a58-4ef6-bfb5-0156650932e7
raw_long_data_comp_2 = DataFrame(
	names = repeat(raw_counts.transcript_id[map_component_back_to_raw],length(2:99)),
	var = stack(raw_counts[map_component_back_to_raw,2:99])[:,1],
	val = stack(raw_counts[map_component_back_to_raw,2:99])[:,2]
);

# ╔═╡ 6112be32-1cff-4eda-8a9c-b2bf350e1e0e
map_processed_back_to_raw = findall(x->x in processed_counts.transcript_id,raw_counts.transcript_id);

# ╔═╡ 67bd2eef-b716-41f0-89b0-4edc07c88948
raw_long_data = DataFrame(
	names = repeat(raw_counts.transcript_id[map_processed_back_to_raw],length(2:99)),
	var = stack(raw_counts[map_processed_back_to_raw,:][:,2:99])[:,1],
	val = stack(raw_counts[map_processed_back_to_raw,:][:,2:99])[:,2]
);

# ╔═╡ e6670e4e-e11a-4603-a4cb-fea6bc603933
begin
	test = zeros(Int,140)
	for i in 1:140 
		test[i] = length(filter(x->x == processed_counts.transcript_id[components[2]][i],raw_counts.transcript_id))
	end
end

# ╔═╡ 72849afd-0f69-4b1b-ba2f-5d7567a6b682
##node that is repeated in network!
prob = processed_counts.transcript_id[components[2]][86]

# ╔═╡ 0c47136e-c65c-421c-bdec-318c61a75282
raw_counts[findall(raw_counts.transcript_id .== prob),:]

# ╔═╡ b46d29a7-6615-4499-8cb0-0a10fb513443
yodel = countmap(raw_counts.transcript_id)

# ╔═╡ a65a8563-ec0e-4627-8ca0-7271a03db258
collect(keys(yodel))[values(yodel).!==1]

# ╔═╡ adc42541-74d1-4ef4-a231-67630241e9d3
yodel["ENST00000335426.4"]

# ╔═╡ Cell order:
# ╠═63f79e1a-b0d8-11ed-1669-a770b3844775
# ╠═5de6c623-1c96-41ba-9bda-534bfbeb8928
# ╠═170965d2-0728-4baf-950e-f2b56739427d
# ╠═02211a8a-87a4-4bd5-acf8-74cd27b902c6
# ╠═5e6beb8d-47fc-4a1b-926e-ae17a51b2ee4
# ╠═db25a7ae-26b0-40da-97d3-06edb70c70d9
# ╠═53770902-0cee-4880-b480-e475e158a369
# ╠═0ed6a565-9abc-428b-bf0f-81be4a0b1f9c
# ╠═392062f3-e17f-4532-9337-692c229944a1
# ╠═96e7e7aa-ec6a-44fd-8914-70c3f44c771f
# ╠═f179209c-cb71-437c-95c6-ba39513ee394
# ╠═df684d90-6453-41f8-b0b2-7ca47a4f9ce3
# ╠═c97ea9f4-c179-4b8e-9abe-d180b75e1abb
# ╠═05cec146-5a7b-4bc6-8797-d7f89f5b8b8b
# ╠═02d67238-533e-4c36-85c4-f828f130e340
# ╠═7882c063-a58f-4747-a3b0-54c4a756bcf4
# ╠═7e828b37-d757-4c23-bf21-aa821e776291
# ╠═69a0ad82-e6e7-4832-b393-95ff5561722e
# ╠═67bd2eef-b716-41f0-89b0-4edc07c88948
# ╠═a854f969-4a58-4ef6-bfb5-0156650932e7
# ╠═d651be50-9d95-49fc-bc70-dfc119947a6c
# ╠═6112be32-1cff-4eda-8a9c-b2bf350e1e0e
# ╠═e6670e4e-e11a-4603-a4cb-fea6bc603933
# ╠═72849afd-0f69-4b1b-ba2f-5d7567a6b682
# ╠═0c47136e-c65c-421c-bdec-318c61a75282
# ╠═a65a8563-ec0e-4627-8ca0-7271a03db258
# ╠═adc42541-74d1-4ef4-a231-67630241e9d3
# ╠═b46d29a7-6615-4499-8cb0-0a10fb513443
