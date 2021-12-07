#need to activate Revise here for now (because Pluto?)
using Revise, JLD2
### Include all source files TODO make this occur more fluently and automatically by creating a package, and using Revise
cwd = ENV["JULIA_PROJECT"]
#define structure for run_parameters (nb... at this stage has to be rerun if values change). Also needs to be defined before including any package that depends on it
struct RunParameters
    test_name::String   
    page_name::String
    website_dir::String
    expression_cutoff::Int
    norm_method::String
    variance_percent::Float64
    coexpression::String
    threshold::Float64
    threshold_method::String
    null_model_size::Int
    func_annotate::Bool
    visualise::Bool
    graphlet_counting::Bool
    wgcna::Bool
end


for src in filter(x->endswith(x,".jl"),readdir("src"))
    if(!occursin("Initialisation.jl",src))
        includet("$cwd/src/"*src)
    end
end

macro name(arg)
    x = string(arg)
    quote
        $x
    end
end

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

function cache_save(file_name::String,save_object::T) where T<:Pair
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

#set up distributed workers
#first clean to make sure there are no stray workers already around
if(length(workers())!=Threads.nthreads())
    rmprocs(workers())
    #add workers equal to the number of available cpus  
    addprocs(Threads.nthreads())
    #addprocs(8)
    #@everywhere include("src/CoexpressionMeasures.jl")
    #@everywhere include("src/CountingFunctions.jl")
    #@everywhere include("src/NullModel.jl")

end
