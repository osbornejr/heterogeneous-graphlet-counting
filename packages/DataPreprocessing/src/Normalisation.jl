using RCall, DataFrames, LinearAlgebra, Statistics, Cairo, Gadfly, Compose ##order is important!!
function library_size_normalisation(raw_counts::Union{DataFrame,Array},method::String)
#library size normalisation methods
#
#TODO make all methods have consistent names, and document these
#These methods all scale across samples to normalise for different library sizes. 
#Methods are all drawn from R packages such as edgeR, limma and EBSeq. 
#Most methods are specifically designed for DE analysis, and may not be important or have adverse affects on GCE analysis.  
#Note, library size methods should only be applied to raw counts, not TPM, FPKM etc as these have already adjusted for libraery size. 
#TPM(transcripts per million) from RSEM is essentially a Total Count library normalisation adjusted to be per million instead of a strict proportion.
    
#fix to fit R package:
    if(method=="upper_quartile")
        method="upperquartile"
    end
    if(method=="total_count")
        #in edgeR this will set norm factors to 1, i.e. allow just to noramlise by straight library size
        method="none"
    end
    #ensure input method to R is valid!
    avail = ["quantile","median","upperquartile","TMM","TMMwsp","none"]
    avail_print = ["quantile","median","upper_quartile","TMM","TMMwsp","total_count"]
    if (!(method in avail) )
        @error "specified normalisation method is not available. Choose from $(avail_print)"
        return Nothing
    end
    
    raw_counts=Array(raw_counts)
    
    #Clear R workspace. May be better implemented as a function?
    R"""
    sapply(names(sessionInfo()$otherPkgs),function(pkg) detach(paste0('package:',pkg),character.only =T,force = T));rm(list=ls())
    """
    @rput raw_counts
    @rput method
    
    R"""
    library(edgeR,quietly=T)
    library(EBSeq,quietly = T)
 
    if(method=="quantile")
    {
      ## Quantile method:
      # simple, rough method that uses the limma package (a dependency of edgeR)
      norm_counts=normalizeQuantiles(raw_counts) 
    }
    
    if(method=="median")
    {
      ## Median method:
      # simple, rough method that uses the EBSeq package
      norm_counts=GetNormalizedMat(raw_counts,MedianNorm(raw_counts,alternative=TRUE))
    }

    if(method %in% c("TMM","TMMwsp","upperquartile","none")) ##uses edgeR?
    {
        list=calcNormFactors(DGEList(raw_counts),method=method)
        ##scale factor to keep numbers within bounds of original data (mean of lib_size)
        scale_factor = mean(list$samples[,2])
        for (i in 1:nrow(list$samples)){list$counts[,i]=scale_factor*list$counts[,i]/(list$samples[i,2]*list$samples[i,3])}
        norm_counts=list$counts
    }
    """
    @rget norm_counts
    return norm_counts
end


# across sample normalisation methods

#pca
"""
    pca(data)
perform PCA using SVD\n
inputs:\n
    - data: M x N matrix of input data. (M dimensions, N trials)\n
outputs:\n
    - PCs: each column is a principle component\n
    - D: N x 1 matrix of variances\n
    - Z: the original data X projected onto the PCs. Z=PCs*X
"""
function pca(data::AbstractArray{T,2}) where T<:Real
   
#Here we compute PCA following the SVD method outlined here: https://medium.com/@jonathan_hui/machine-learning-singular-value-decomposition-svd-principal-component-analysis-pca-1d45e885e491
#Input should be an M x N array that features M observations over N trials.
# Slightly modifies and explains function given by comment at https://discourse.julialang.org/t/how-to-get-the-principal-components-and-variances-from-multivariatestats/15843/2
    #remove rows with zero variance
    v_data=data[vec(std(data,dims=2).!=0),:]
    #standardise data
    X = (v_data .- mean(v_data, dims=2)) #./ std(v_data,dims=2)
    #substitute so that Y'*Y is cov(X)
    Y = X' ./ sqrt(T(size(X,2)-1))
    #Compute (truncated) SVD of Y, which gives PCs of X (eigenvectors of cov(X)). See Shlens (2014) tutorial.
    U,S,V = svd(Y)
        Σ = diagm(0=>S)
    #by same token, the square of the singular values diagonal matrix gives us the eigenvalues (variances) of cov(X)
    v = Σ * Σ
        #sort variance values and PCs 
        indexList = sortperm(diag(v); rev=true)
    PCs = map(x->V[:,x], indexList)
    #flatten components into one array (desirable? maybe label?)
    PCs = hcat(PCs...)
    #use method outlined at https://blog.bioturing.com/2018/06/14/principal-component-analysis-explained-simply/ 
    #to generate a per-sample value in each component 
    #(experimental, need a reference/verification plus testing. Could be nice to name component columns) 
    #per_sample=X'*PCs
    #
    #Finally, provide the original data projected onto the PCs 
    Z = V'*v_data

#Note that following the SVD method, we just need to take the transpose of V to get the PCs. This gives the correct (expected) dimensions for each component rather than the error from before. 
        return V', S.^2, Z
end



##NOTE as of recently, pca above has been tidied up, so now these plot functions expect input to be Z, the projection of the data onto the PCs, with each column representing a PC.
function pca_plot(Z::AbstractArray{T,2},grid_dims::Int) where T<:Real
#this function allows visualisation of the PCA decomposition. Grid dimensions sets the number of principal components that will be compared. A grid_dim of 1 implies a simple biplot of the first 2 principal  components.
    data=Z#DataFrame(PCs,:auto)
    #store plots in array
    plots=Array{Plot}(undef,grid_dims,grid_dims)
    for i in 1:grid_dims
        for j in 1:grid_dims
            plots[i,j]=Gadfly.plot(data,x=data[j+i,:],y=data[i,:],Guide.xlabel(string("PC ",j+i)),Guide.ylabel(string("PC ",i)),Theme(key_position=:none)) 
        end
    end
    Gadfly.set_default_plot_size(20cm, 20cm)
    Grid=gridstack(plots)
    return Grid
end 

function pca_plot(Z::AbstractArray{T,2},grid_dims::Int,filename::String) where T<:Real
#this function allows visualisation of the PCA decomposition. Grid dimensions sets the number of principal components that will be compared. A grid_dim of 1 implies a simple biplot of the first 2 principal  components.
#This version is to save plot to file, at given path filename
    data=Z#DataFrame(PCs,:auto)
    #store plots in array
    plots=Array{Plot}(undef,grid_dims,grid_dims)
    for i in 1:grid_dims
        for j in 1:grid_dims
            plots[i,j]=Gadfly.plot(data,x=data[j+i,:],y=data[i,:],Guide.xlabel(string("PC ",j+i)),Guide.ylabel(string("PC ",i)),Theme(key_position=:none)) 
        end
    end
    Gadfly.set_default_plot_size(20cm, 20cm)
    Grid=gridstack(plots)
    draw(SVG(filename),Grid)
end 

# transformations (log2, VST, etc)

# sample box plot

function boxplot(dataframe::DataFrame,filename::String,groupdict::Dict{String,String})
    longform = stack(dataframe,variable_name = "sample")
    longform = longform[longform[!,:value].!=0,:]
    insertcols!(longform,"log_value"=>log2.(longform[!,:value]))
    #Add column for color group
    insertcols!(longform,"group"=>[groupdict[samp] for samp in longform[!,:sample]])
    p = plot(longform, x = "sample", y = "log_value", Geom.boxplot(suppress_outliers = true),Guide.xticks(label=false),color=:group);
    draw(SVG(filename),p)
end

#Another boxplot method that defaults to colouring based on sample
function boxplot(dataframe::DataFrame,filename::String)
    longform = stack(dataframe,variable_name = "sample")
    longform = longform[longform[!,:value].!=0,:]
    insertcols!(longform,"log_value"=>log2.(longform[!,:value]))
    #Add column for color group
    p = plot(longform, x = "sample", y = "log_value", Geom.boxplot(suppress_outliers = true),Guide.xticks(label=false),color=:sample);
    draw(SVG(filename),p)
end

#Another boxplot method that defaults to colouring based on sample and outputs plot object
function boxplot(dataframe::DataFrame)
    longform = stack(dataframe,variable_name = "sample")
    #Add column for color group
    p = plot(longform, x = "sample", y = "value", Geom.boxplot(suppress_outliers = true),Guide.xticks(label=false),Theme(key_position= :none,grid_line_width=0mm),color=:sample);
    return p
end

function histogram(data::Vector,filepath;title::String="",xaxis::String="",yaxis::String="frequency")
    p = plot(x = log2.(data), Geom.histogram(bincount = 100,density = false),Guide.xlabel(xaxis),Guide.ylabel(yaxis),Guide.title(title));
    draw(SVG(filepath),p)
end
function histogram(data::DataFrame,plotcol::Symbol,colourcol::Symbol,filepath;title::String="",xaxis::String="",yaxis::String="frequency")
    p = plot(data,x = plotcol, Geom.histogram(bincount = 100,density = false),Guide.xlabel(xaxis),Guide.ylabel(yaxis),Guide.title(title),color=colourcol);
    draw(SVG(filepath),p)
end

function histogram(data::DataFrame,plotcol::Symbol,colourcol::Symbol;title::String="",xaxis::String="",yaxis::String="frequency")
    p = plot(data,x = plotcol, Geom.histogram(bincount = 100,density = false),Guide.xlabel(xaxis),Guide.ylabel(yaxis),Guide.title(title),color=colourcol);
    return p
end

  
