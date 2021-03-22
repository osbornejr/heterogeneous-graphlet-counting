using Statistics, StatsBase, DataFrames, InformationMeasures, LinearAlgebra, RCall,Distributed,ProgressMeter

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
	if (method=="PID")
		return partial_information_decomposition(data; discretizer = "uniform_width", estimator = "maximum_likelihood", mi_base = 2)
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

#aggregator function
function t2(d1,d2)
	append!(d1,d2)
	d1
		
end
function get_unique_values(data,i,j)		#set unique vectors
	ux_vector = Vector{Float64}(undef,size(data,1))
	##calculate unique scores for every triplet involving i and j 
	for k in 1:size(data,1)
		ux_vector[k] = get_partial_information_decomposition(data[j,:],data[k,:],data[i,:])["unique_1"]
	end
	return sum(ux_vector)
end


function partial_information_decomposition(data; discretizer = "uniform_width", estimator = "maximum_likelihood", mi_base = 2,distributed::Bool = false )
	## set up bins in advance, we can then calculate probabilities within triple loop 
	nvars, nvals = size(data)
	matrix = zeros(nvars,nvars)
	if (distributed = true)
		##TODO need shared array for larger data matrices?
		uniques = @showprogress @distributed (t2) for t in [(x,y) for x in 1:nvars, y in 1 : nvars]
			[get_unique_values(data,first(t),last(t))]
		end
		uniques = reshape(uniques,nvars,nvars)

	else
		for i in 1 : nvars, j in i : nvars
			@info "Calculating for gene pair $i, $j..."
			#reset unique vectors
			ux_vector = Vector{Float64}(undef,nvars)
			uy_vector = Vector{Float64}(undef,nvars)
			##calculate unique scores for every triplet involving i and j 
			for k in 1: nvars
				ux_vector[k] = get_partial_information_decomposition(data[j,:],data[k,:],data[i,:])["unique_1"]
				uy_vector[k] = get_partial_information_decomposition(data[i,:],data[k,:],data[j,:])["unique_1"]
			end
		end
	end
	MI = mutual_information(data)
	#add unique_xy and unique_yx together, maintaining symmetry in matrix
	PUC = UpperTriangular(uniques)'+LowerTriangular(uniques)+UpperTriangular(uniques)+LowerTriangular(uniques)'
	# normalise by MI values
	PUC = PUC./MI
	return PUC
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


	
