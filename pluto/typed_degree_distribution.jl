### A Pluto.jl notebook ###
# v0.19.14

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
	#using DataPreprocessing
	using NetworkConstruction
	using GraphletCounting
	#using GraphletAnalysis
	using ProjectFunctions
end

# ╔═╡ 5de6c623-1c96-41ba-9bda-534bfbeb8928
# ╠═╡ show_logs = false
load_config(cwd*"/config/run-files/GSE68559.yaml")

# ╔═╡ 170965d2-0728-4baf-950e-f2b56739427d
# ╠═╡ show_logs = false
adj_matrix,network_counts,vertexlist,edgelist = network_construction();

# ╔═╡ 02211a8a-87a4-4bd5-acf8-74cd27b902c6
dd = GraphletCounting.typed_degree_distribution(vertexlist,edgelist);

# ╔═╡ 5e6beb8d-47fc-4a1b-926e-ae17a51b2ee4
data = DataFrame(coding = map(x->DefaultDict(0,x)["coding"],dd),noncoding = map(x->DefaultDict(0,x)["noncoding"],dd))

# ╔═╡ bd3d754d-9687-4ac2-bccc-b16fca7dfcf5
plot(data[vertexlist.=="noncoding",:], x = "noncoding",Geom.histogram,Guide.title("Typed degree distribution"),style(point_size = 0.4mm))

# ╔═╡ 732e425b-6549-48c5-bf8c-e196d977356b
plot(data[vertexlist.=="noncoding",:], x = "noncoding",Geom.histogram,Guide.title("Typed degree distribution"),style(point_size = 0.4mm))

# ╔═╡ 0502e056-3059-4941-84d0-fe7da897bc48
data[vertexlist.=="coding",:]

# ╔═╡ Cell order:
# ╠═63f79e1a-b0d8-11ed-1669-a770b3844775
# ╠═5de6c623-1c96-41ba-9bda-534bfbeb8928
# ╠═170965d2-0728-4baf-950e-f2b56739427d
# ╠═02211a8a-87a4-4bd5-acf8-74cd27b902c6
# ╠═5e6beb8d-47fc-4a1b-926e-ae17a51b2ee4
# ╠═bd3d754d-9687-4ac2-bccc-b16fca7dfcf5
# ╠═732e425b-6549-48c5-bf8c-e196d977356b
# ╠═0502e056-3059-4941-84d0-fe7da897bc48
