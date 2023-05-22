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
    #ensure input method to R is valid!
    avail = ["quantile","median","upperquartile","TMM","TMMswp","RLE"]
    if (!(method in avail) )
        @error "specified normalisation method is not available. Choose from $(avail)"
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

    if(method=="DESeq2")
    {  
    ##DESeq2 method (NOT WORKING; is this necessary in any case, considering assumptions of DESeq2 methods?):
     #compute geometric mean for each gene across samples. Doesn't handle zero values well (unclear how they are handled in package)
      geo_mean_total_homo_sap <- total_homo_sap %>%
        data_matrix() %>% 
        .^(1/length(.)) %>% 
        rowProds() %>%
        as_tibble()
    }
    
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
      norm_counts=GetNormalizedMat(raw_counts,MedianNorm(raw_counts))
    }

    if(method %in% c("TMM","TMMwsp","upperquartile","RLE")) ##uses edgeR?
    {
        list=calcNormFactors(DGEList(raw_counts),method=method)
        for (i in 1:nrow(list$samples)){list$counts[,i]=1e6*list$counts[,i]/(list$samples[i,2]*list$samples[i,3])}
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
    - PC: each column is a principle component\n
    - V: M x 1 matrix of variances\n
"""
function pca(data::AbstractArray{T,2}) where T<:Real
   
#Here we compute PCA following the SVD method outlined here: https://medium.com/@jonathan_hui/machine-learning-singular-value-decomposition-svd-principal-component-analysis-pca-1d45e885e491
#Input should be an M x N array that features M observations over N trials.
# Slightly modifies and explains function given by comment at https://discourse.julialang.org/t/how-to-get-the-principal-components-and-variances-from-multivariatestats/15843/2
    #remove rows with zero variance
    v_data=data[vec(std(data,dims=2).!=0),:]
    #standardise data
    X = (v_data .- mean(v_data, dims=2)) ./ std(v_data,dims=2)
    #substitute so that Y'*Y is cov(X)
    Y = X' ./ sqrt(T(size(X,2)-1))
    #Compute SVD of Y, which gives PCs of X (eigenvectors of cov(X)). See Shlens (2014) tutorial.
    U,S,PC = svd(Y)
        Σ = diagm(0=>S)
    #by same token, the square of the singular values diagonal matrix gives us the eigenvalues (variances) of cov(X)
    V = Σ * Σ
        #sort variance values and PCs 
        indexList = sortperm(diag(V); rev=true)
    PCs = map(x->PC[:,x], indexList)
    #flatten components into one array (desirable? maybe label?)
    PCs = hcat(PCs...)
    #use method outlined at https://blog.bioturing.com/2018/06/14/principal-component-analysis-explained-simply/ 
    #to generate a per-sample value in each component 
    #(experimental, need a reference/verification plus testing. Could be nice to name component columns) 
    per_sample=X'*PCs
        return PCs, diag(V)[indexList], per_sample
end


function pca_plot(PCs::AbstractArray{T,2},grid_dims::Int) where T<:Real
#this function allows visualisation of the PCA decomposition. Grid dimensions sets the number of principal components that will be compared. A grid_dim of 1 implies a simple biplot of the first 2 principal  components.
    data=DataFrame(PCs,:auto)
    #store plots in array
    plots=Array{Plot}(undef,grid_dims,grid_dims)
    for i in 1:grid_dims
        for j in 1:grid_dims
            plots[i,j]=Gadfly.plot(data,x=data[:,j+i],y=data[:,i],Guide.xlabel(string("PC ",j+i)),Guide.ylabel(string("PC ",i)),Theme(key_position=:none)) 
        end
    end
    Gadfly.set_default_plot_size(20cm, 20cm)
    Grid=gridstack(plots)
    return Grid
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
    longform = longform[longform[!,:value].!=0,:]
    insertcols!(longform,"log_value"=>log2.(longform[!,:value]))
    longform = longform[longform[!,:log_value].>-10,:]
    longform = longform[longform[!,:log_value].<10,:]
    #Add column for color group
    p = plot(longform, x = "sample", y = "log_value", Geom.boxplot(suppress_outliers = true),Guide.xticks(label=false),Theme(key_position= :none),color=:sample);
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
