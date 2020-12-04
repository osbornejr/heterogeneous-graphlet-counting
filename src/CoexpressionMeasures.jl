using Statistics, StatsBase, DataFrames, InformationMeasures, LinearAlgebra, RCall

function coexpression_measure(data::Union{AbstractDataFrame,AbstractArray},method::String)
	if (method =="pearson")
		return cor(data')
	end

	if (method=="spearman")
		return corspearman(data')
	end

	if (method== "kendall")
		return corkendall(data')
	end

	if (method=="mutual_information")
		return mutual_information(data; discretizer = "uniform_width", estimator = "maximum_likelihood", mi_base = 2)
	end
	if (method=="pcit")
		return pcit(data)
	end
end


#modified form of function described on InformationMeasures.jl github page. Faster than old method because bin calculation isn't repeated for each variable for every element. Main changes are orientation (variables as rows) and outputting as a symmetric matrix rather than a one dimensional array.
function mutual_information(data; discretizer = "uniform_width", estimator = "maximum_likelihood", mi_base = 2)

	nvars, nvals = size(data)

	bin_ids = zeros(Int, (nvars, nvals))
	nbins = Int(round(sqrt(nvals)))
	#mis = zeros(binomial(nvars, 2))
	matrix = zeros(nvars,nvars)
	for i in 1 : nvars
		get_bin_ids!(view(data,i,1:nvals), discretizer, nbins, view(bin_ids, i, 1:nvals))
	end

	#index = 1
	for i in 1 : nvars, j in i : nvars
		f = get_frequencies_from_bin_ids(view(bin_ids,i,1:nvals), view(bin_ids,j,1:nvals), nbins, nbins)
		p = get_probabilities(estimator, f) 
		#mis[index] = apply_mutual_information_formula(p, sum(p, dims = 1), sum(p, dims = 2), mi_base)
		matrix[i,j] = apply_mutual_information_formula(p, sum(p, dims = 1), sum(p, dims = 2), mi_base)

		#index += 1
	end
	#copy upper triangle to lower triangle
	matrix = matrix + matrix'
	matrix[diagind(matrix)]=matrix[diagind(matrix)]./2
	return matrix

end

function pcit(data::AbstractArray)
	#implementing the R package PCIT here for now. In the long run it might be better/faster/more desirable to implement our own version in Julia?
	# outputs a correlation matrix with those values deemed insignificant by the pcit algorithm set to 0.
	@rput data
	R"""
	library(PCIT)
	#just using pearson cor as basis for now (could use nonlinear correlation method here? unclear from PCIT documentation)
	c = cor(t(data))
	results = pcit(c)
	sig = idx(results)
	nonsig = idxInvert(nrow(c),sig)
	#at the moment, just choose significant interactions in binary (i.e. only for unweighted networks)
	c[nonsig] = 0
	"""
	@rget c
	return c
end
function consensus_measure(data::AbstractArray;methods::AbstractVector{String})
	nvars, nvals = size(data)
	consensus_matrix = ones(nvars,nvars)
	threshold = 0.95
	for method in methods
		similarity_matrix = coexpression_measure(data,method)
		adjacency_matrix = adjacency(similarity_matrix,threshold)
		consensus_matrix = consensus_matrix.*adjacency_matrix
	end
	return BitArray(consensus_matrix)
end


	
