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
	Guide.ylabel("frequency"))

# ╔═╡ 0ed6a565-9abc-428b-bf0f-81be4a0b1f9c
components = NetworkConstruction.network_components(adj_matrix) 

# ╔═╡ 392062f3-e17f-4532-9337-692c229944a1
 processed_data = DataPreprocessing.data_from_dataframe(processed_counts,"data")

# ╔═╡ 96e7e7aa-ec6a-44fd-8914-70c3f44c771f
dodgy_component_data = processed_data[components[2],:]

# ╔═╡ 7e828b37-d757-4c23-bf21-aa821e776291
processed_long_data = DataFrame(
	names = repeat(processed_counts.transcript_id,length(2:99)),
	var = stack(processed_counts[:,2:99])[:,1],
	val = stack(processed_counts[:,2:99])[:,2]
);

# ╔═╡ 69a0ad82-e6e7-4832-b393-95ff5561722e
processed_long_data_comp_2 = DataFrame(
	names = repeat(processed_counts.transcript_id[components[2]],length(2:99)),
	var = stack(processed_counts[components[2],2:99])[:,1],
	val = stack(processed_counts[components[2],2:99])[:,2]
);

# ╔═╡ a854f969-4a58-4ef6-bfb5-0156650932e7
raw_long_data_comp_2 = DataFrame(
	names = repeat(raw_counts.transcript_id[components[2]],length(2:99)),
	var = stack(raw_counts[components[2],2:99])[:,1],
	val = stack(raw_counts[components[2],2:99])[:,2]
);

# ╔═╡ f179209c-cb71-437c-95c6-ba39513ee394
Gadfly.plot(raw_long_data_comp_2,x=:var,y=:val,group=:names,Geom.line,Guide.xticks(label=false),Theme(key_position = :none))


# ╔═╡ d651be50-9d95-49fc-bc70-dfc119947a6c
map_component_back_to_raw = findall(x->x in processed_counts.transcript_id[components[2]],raw_counts.transcript_id)

# ╔═╡ e6670e4e-e11a-4603-a4cb-fea6bc603933
begin
	test = zeros(Int,140)
	for i in 1:140 
		test[i] = length(filter(x->x == processed_counts.transcript_id[components[2]][i],raw_counts.transcript_id))
	end
end

# ╔═╡ 72849afd-0f69-4b1b-ba2f-5d7567a6b682
prob = processed_counts.transcript_id[components[2]][86]

# ╔═╡ 0c47136e-c65c-421c-bdec-318c61a75282
raw_counts[findall(raw_counts.transcript_id .== prob),:]

# ╔═╡ a65a8563-ec0e-4627-8ca0-7271a03db258


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
# ╠═7e828b37-d757-4c23-bf21-aa821e776291
# ╠═69a0ad82-e6e7-4832-b393-95ff5561722e
# ╠═a854f969-4a58-4ef6-bfb5-0156650932e7
# ╠═d651be50-9d95-49fc-bc70-dfc119947a6c
# ╠═e6670e4e-e11a-4603-a4cb-fea6bc603933
# ╠═72849afd-0f69-4b1b-ba2f-5d7567a6b682
# ╠═0c47136e-c65c-421c-bdec-318c61a75282
# ╠═a65a8563-ec0e-4627-8ca0-7271a03db258
