module ProjectFunctions
export cwd,cache_save, cache_load,  @name
using DataPreprocessing, NetworkConstruction,GraphletCounting ,GraphletAnalysis

include("RunFunctions.jl")

cwd = ENV["PWD"]
#define structure for run_parameters (nb... at this stage has to be rerun if values change). Also needs to be defined before including any package that depends on it
#struct RunParameters
#    test_name::String   
#    page_name::String
#    website_dir::String
#    expression_cutoff::Int
#    norm_method::String
#    variance_percent::Float64
#    coexpression::String
#    threshold::Float64
#    threshold_method::String
#    null_model_size::Int
#    func_annotate::Bool
#    visualise::Bool
#    graphlet_counting::Bool
#    wgcna::Bool
#end

function distributed_setup(inclusions::Vector{Symbol})
    distributed_setup(inclusions...)
end

function distributed_setup(inclusions::Symbol...)

    ##set up for distributed mode

    #first clean to make sure there are no stray workers already around
    rmprocs(workers())
    #add workers equal to the number of available cpus      
    addprocs(Threads.nthreads();exeflags="--project=$cwd")
    #addprocs(8)
    #@everywhere inclusions
    for inc in inclusions
        eval(macroexpand(Distributed,quote @everywhere using $(inc) end))
    end
end

function distributed_code_load(inclusions::Symbol...)
    for package in inclusions
        eval(macroexpand(Distributed,quote @everywhere using $(package) end))
    end
end

macro name(arg)
    x = string(arg)
    quote
        $x
    end
end

#aggregator function
function t2(d1,d2)
    append!(d1,d2)
    d1
end

##allows a function to be broadcast and splat each argument
splat(f) = args->f(args...)

function cache_save(file_name::String,save_objects::Array{<:Pair,1})
    ##shorthand way to use JLD2 method to save files to disk
    if (split(file_name,".")[end]!="jld2")
        throw(DomainError(file_name,"output file must be a valid .jld2 file"))
    end
    JLD2.jldopen(file_name,"w") do file
        for o in save_objects
            file[first(o)] = last(o)
        end
    end;
end

function cache_save(file_name::String,save_object::T) where T<:Pair{String,}
    ##shorthand way to use JLD2 method to save files to disk
    if (split(file_name,".")[end]!="jld2")
        throw(DomainError(file_name,"output file must be a valid .jld2 file"))
    end
    JLD2.jldopen(file_name,"w") do file
            file[first(save_object)] = last(save_object)
    end;
end

function cache_load(file_name::String,load_object::String)
    #load specific object from jld2 file
    if (split(file_name,".")[end]!="jld2")
        throw(DomainError(file_name,"input file must be a valid .jld2 file"))
    end
    JLD2.jldopen(file_name) do file
        file[load_object]
    end
end


    
end # module
