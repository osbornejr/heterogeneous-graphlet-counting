using LinearAlgebra,RCall
function adjacency(data::AbstractArray,threshold::Float64)
	sim_matrix = copy(data)
	sim_matrix[diagind(sim_matrix)].= 0
	sim_matrix[broadcast(abs,sim_matrix).<threshold].=0
	sim_matrix[broadcast(abs,sim_matrix).>threshold].=1
	sim_matrix=BitArray(sim_matrix)
	return sim_matrix
end
function empirical_dist_zero_adjacency(sim_matrix::AbstractArray,prob::Float64)
	##run each row of similarity matrix through its own empirical cdf to determine probability 
	dists = hcat(map.(ecdf.(eachrow(sim_matrix)),eachrow(sim_matrix))...)
	##alternative method to ignore 0 values in ecdfs...
	#dists = hcat(map.(ecdf.(filter.(x->x!=0,eachrow(sim_matrix))),eachrow(sim_matrix))...)
	##find which are significant for each row and column, then combine
	adj = (dists.>prob).*(dists'.>prob)
	adj[diagind(adj)].= 0
	return adj
end
function empirical_dist_adjacency(sim_matrix::AbstractArray,prob::Float64)
	##run each row of similarity matrix through its own empirical cdf to determine probability 
	#dists = hcat(map.(ecdf.(eachrow(sim_matrix)),eachrow(sim_matrix))...)
	##alternative method to ignore 0 values in ecdfs...
	dists = hcat(map.(ecdf.(filter.(x->x!=0,eachrow(sim_matrix))),eachrow(sim_matrix))...)
	##find which are significant for each row and column, then combine
	adj = (dists.>prob).*(dists'.>prob)
	adj[diagind(adj)].= 0
	return adj
end

function top_adjacency(sim_matrix::AbstractArray,top_scores::Int)
	sim_matrix[diagind(sim_matrix)].= 0
	##sort each row of similarity matrix
	sorted = map(sort,eachrow(broadcast(abs,sim_matrix)))
	## get lowest included score out of top scores for each row
	lims = getindex.(sorted,Ref(size(sim_matrix,1)-top_scores))
	cuts = broadcast(abs,sim_matrix).>lims
	##check against transpose to see which entries agree for BOTH transcripts
	adj = cuts.*cuts'
	return adj
end

function edgelist_from_adj(adjacency_matrix::AbstractArray)
	edgelist=Array{Pair}(undef,sum(UpperTriangular(adjacency_matrix)))
	count=0
	for (i,row) in enumerate(eachrow(UpperTriangular(adjacency_matrix)))
       		for j in 1:size(row,1)
			if (row[j]==1)
				count=count+1
				edgelist[count]=Pair(i,j)
       			end
       		end
       end
	return edgelist
end


function wgcna_network(data::AbstractArray,transcript_types::Array{String})
	#transpose for WGCNA
	data = data'
	@rput data	
	@rput transcript_types 
	R"""
	library(WGCNA)	
	library(dendextend)
	library(ggdendro)
	library(tidyverse)
	##save as one pdf
	pdf('wgcna.pdf')
	# Choose a set of soft-thresholding powers 
	powers = c(c(1:10), seq(from = 12, to=20, by=2))
	# Call the network topology analysis function 
	sft = pickSoftThreshold(data, powerVector = powers, verbose = 5) 
	# Scale-free topology fit index as a function of the soft-thresholding power 
	#pdf('powers.pdf')
	plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence")) 
	text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,col="red") 
	# this line corresponds to using an R^2 cut-off of h 
	abline(h=0.90,col="red") 
	# Mean connectivity as a function of the soft-thresholding power 
	plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power) for ICC283",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity")) 
	text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red") 
	#dev.off()
	softPower = 9; 
	adjacency = adjacency(data, power = softPower) 
	TOM=TOMsimilarity(adjacency) 
	dissTOM=1-TOM 
	
	##auto-option
	#bwnet = blockwiseModules(data, maxBlockSize = 4000, power = 8, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, saveTOMs = TRUE, saveTOMFileBase = "data-TOM-blockwise", verbose = 3) 
	
	## Dendro
	geneTree = hclust(as.dist(dissTOM), method = "average");  
	#dend object (for dendextend etc)
	genedend = hang.dendrogram(as.dendrogram(geneTree),hang = 0.04)
	#set colours of leaves
	#genedend<- branches_attr_by_labels(genedend,(1:length(transcript_types))[transcript_types=="coding"],attr = c("col"), type = c("all")) 
	label_colours = ifelse(labels(genedend)%in%(1:length(transcript_types))[transcript_types=="coding"],"red","blue")
	genedend<- assign_values_to_leaves_edgePar(genedend, value = label_colours, edgePar = "col") 
	# Plot the resulting clustering tree (dendrogram)  
	#pdf('dendro.pdf')
	#plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04);  
	plot(dendrapply(genedend,noLabel))
	#p<- plot_ggdendro(hcdata, direction="tb",expand.y = 0.2)
	#dev.off()

	## Modules
	minModuleSize = 30;  
	# Module identification using dynamic tree cut:  
	dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);  
	table(dynamicMods)  
	# Convert numeric lables into colors  
	dynamicColors = labels2colors(dynamicMods)  
	table(dynamicColors)  
	# Plot the dendrogram and colors underneath  
	#pdf('colours.pdf')
	plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",  dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")  
	#dev.off()
	
	# Calculate eigengenes  
	MEList = moduleEigengenes(data, colors = dynamicColors)  
	MEs = MEList$eigengenes  
	# Calculate dissimilarity of module eigengenes  
	MEDiss = 1-cor(MEs);  
	# Cluster module eigengenes  
	METree = hclust(as.dist(MEDiss), method = "average");  
	# Plot the result  
	#pdf("eigengenes.pdf")
	plot(METree, main = "Clustering of module eigengenes",  xlab = "", sub = "")  
	MEDissThres = 0.25  
	# Plot the cut line into the dendrogram  
	abline(h=MEDissThres, col = "red")  
	#dev.off()
	# Call an automatic merging function  
	merge = mergeCloseModules(data, dynamicColors, cutHeight = MEDissThres, verbose = 3)  
	# The merged module colors  
	mergedColors = merge$colors;  
	# Eigengenes of the new merged modules:  
	mergedMEs = merge$newMEs;  
	#pdf(file = "eigen-merged-dendro.pdf", wi = 9, he = 6)  
	plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),  c("Dynamic Tree Cut", "Merged dynamic"),  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
	#dev.off()  
	str(mergedMEs)  
	
	
	##Heat map
	# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap 
	plotTOM = dissTOM^7; 
	# Set diagonal to NA for a nicer plot 
	diag(plotTOM) = NA; 
	# Call the plot function 
	#pdf('heatmap.pdf') 
	TOMplot(plotTOM, geneTree, dynamicColors, main = "Network heatmap plot, all genes") 
	dev.off()	 
	 
	"""

end
