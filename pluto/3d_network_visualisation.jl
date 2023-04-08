### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ ca75971e-d54a-11ed-20fe-c374d41bc379
# ╠═╡ show_logs = false
begin
	cwd = ENV["PWD"];	
	import Pkg
	dev = Pkg.develop
	dev(path=cwd*"/packages/DataPreprocessing")
	dev(path=cwd*"/packages/NetworkConstruction")
	dev(path=cwd*"/packages/GraphletCounting")
	dev(path=cwd*"/packages/GraphletAnalysis")
	dev(path=cwd*"/packages/ProjectFunctions")
	using Revise	
	using CommonMark
	using PlutoUI
	using LinearAlgebra
	## load dev packages last
	using DataPreprocessing
	using NetworkConstruction
	using GraphletCounting
	using GraphletAnalysis
	using ProjectFunctions

	using WGLMakie
	using Graphs
	using GraphMakie
	using NetworkLayout
	using JSServe
	Page()
end

# ╔═╡ 7c2b062c-cfa5-4f88-91cf-c776f91f1d44
load_config(cwd*"/config/run-files/GSE68559.yaml") 

# ╔═╡ af6b4890-0db9-4a4d-9ff2-c9fdfa0a72cc
# ╠═╡ show_logs = false
raw_counts,processed_counts,similarity_matrix,adj_matrix,network_counts,vertexlist,edgelist = get_output_data();

# ╔═╡ dbb99ec2-4774-4993-a4b8-dda027dda19f
g = Graph(adj_matrix)

# ╔═╡ 86ca01a5-ffd3-4b91-838d-01c465346593
vertex_colors = replace(vertexlist,"noncoding"=>:red,"coding"=>:blue);

# ╔═╡ a007f64f-e40d-4db3-8e02-25502ea41c51
# ╠═╡ show_logs = false
begin
	set_theme!(backgroundcolor="#121212")
	fig,scene,p = graphplot(g;
		layout=Spring(dim=3,C=1.0),
		node_color = vertex_colors,
		node_size = 5,
		edge_color = :grey,
		edge_width = 0.1,
		figure = (resolution = (1500, 800),)
							
	)
	scene.show_axis =false
	fig
end

# ╔═╡ Cell order:
# ╠═ca75971e-d54a-11ed-20fe-c374d41bc379
# ╠═7c2b062c-cfa5-4f88-91cf-c776f91f1d44
# ╠═af6b4890-0db9-4a4d-9ff2-c9fdfa0a72cc
# ╠═dbb99ec2-4774-4993-a4b8-dda027dda19f
# ╠═a007f64f-e40d-4db3-8e02-25502ea41c51
# ╠═86ca01a5-ffd3-4b91-838d-01c465346593
