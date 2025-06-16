### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ b4cf09da-f855-11ed-0d6e-2b37d5b5d574
# ╠═╡ show_logs = false
begin
	###NOTE: This is a reconstruction of an old weave notebook from late 2020.
cwd = ENV["PWD"];	
	import Pkg
	#Pkg.activate(mktempdir())
    #Pkg.add([
    #    Pkg.PackageSpec(name="DataFrames",version'"1.3"),
	#	Pkg.PackageSpec(name="Revise"),
	#	Pkg.PackageSpec(name="CommonMark"),
	#	Pkg.PackageSpec(name="PlutoUI"),
	#Pkg.PackageSpec(name="DataStructures"),
	#	Pkg.PackageSpec(name="StatsBase"),
	#	Pkg.PackageSpec(name="Gadfly"),
	#	Pkg.PackageSpec(name="LinearAlgebra")
    #])
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

# ╔═╡ 1bf1d48f-d5d1-4a29-955e-e49a02777f4a
md"""
# Examining data variance
## Joel Robertson
## 9th October 2020
"""

# ╔═╡ 1afd118b-12ba-4f9f-b567-2754dcc142c7
raw_counts=DataPreprocessing.read_count_data("$cwd/data/mayank-de-novo/isoforms",method="expected_count");

# ╔═╡ Cell order:
# ╠═b4cf09da-f855-11ed-0d6e-2b37d5b5d574
# ╠═1bf1d48f-d5d1-4a29-955e-e49a02777f4a
# ╠═1afd118b-12ba-4f9f-b567-2754dcc142c7
