module ProjectFunctions
using DataPreprocessing, NetworkConstruction,GraphletCounting ,GraphletAnalysis

include("RunFunctions.jl")
include("CacheFunctions.jl")
include("FetchFunctions.jl")
include("WebFunctions.jl")

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


function cleaner()
    ##remind GC to do its job
    Threads.@threads for i in 1:10000
        zeros(10000)
    end 
    ##notify GC of its duties (including every worker)
    @everywhere GC.gc()
    ## bonus clean for sticky memory that can accumulate (again, for the worker)
    @everywhere ccall(:malloc_trim,Cvoid,(Cint,),0)
end
export cleaner

function distributed_cleanup()
    cleaner()
    rmprocs(workers())
    cleaner()
end
export distributed_cleanup

function distributed_setup(inclusions::Vector{Symbol})
    distributed_setup(inclusions...)
end

function distributed_setup(inclusions::Symbol...;nthreads::Int=0)

    ##set up for distributed mode

    #first clean to make sure there are no stray workers already around
    rmprocs(workers())
    #add workers equal to the number of available cpus if no threads is provided
    if nthreads == 0
        nthreads = Threads.nthreads();
    end

    addprocs(nthreads;exeflags="--project=$cwd")
    #@everywhere inclusions
    for inc in inclusions
        eval(macroexpand(Distributed,quote 
                                        @everywhere using $(inc) 
                                     end))
    end
end
export distributed_setup

function distributed_code_load(inclusions::Symbol...)
    for package in inclusions
        eval(macroexpand(Distributed,quote     
                                        @everywhere using $(package) 
                                    end))
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
export cache_save

function cache_load(file_name::String,load_object::String)
    #load specific object from jld2 file
    if (split(file_name,".")[end]!="jld2")
        throw(DomainError(file_name,"input file must be a valid .jld2 file"))
    end
    @info "Loading $(load_object) data from $(file_name)"
    JLD2.jldopen(file_name) do file
        file[load_object]
    end
end
export cache_load

# macros to load all functions of module into REPL (for development): taken from https://discourse.julialang.org/t/exportall/4970/18
#
#
## Copyright 2021 Clemens Cords
## Created on 26.12.2021 by clem (mail@clemens-cords.com)
##

"""
export all non-temporary, non-imported values (values not having a '#' at the
start of it's symbol when listed via Base.names) in a module

Example:
module MyModule

module SubModule
_sub_variable
end

_variable
end

@make_public MyModule

# exports:
#   MyModule
#   MyModule.SubModule
#   MyModule._variable

See also: @make_public_rec
"""
macro make_public(module_name::Symbol)

    eval(Meta.parse("export " * string(module_name)))

    as_module = eval(module_name)
    @assert as_module isa Module

    for name in names(as_module; all = true)
        if (string(name)[1] != '#')
            #println("export " * string(name))
            as_module.eval(Meta.parse("export " * string(name)))
        end
    end

    return nothing
end
export @make_public

"""
export all non-temporary, non-imported values (values not having a '#' at the
start of it's symbol when listed via Base.names) in a module and all
such values in any submodule, recursively

Example:
module MyModule

module SubModule
_sub_variable
end

_variable
end

@make_public_rec MyModule

# exports:
#   MyModule
#   MyModule.SubModule
#   MyModule.SubModule._sub_variable
#   MyModule._variable

See also: @make_public
"""
macro make_public_rec(module_name::Symbol)

    function make_public_aux(child_name::Symbol, super::Module) ::Nothing

        if (string(child_name)[1] != '#')
            #println("export " * string(child_name))
            super.eval(Meta.parse("export " * string(child_name)))
        end

        child = super.eval(child_name)
        if (child isa Module && child != super)
            for name in names(child; all = true)
                make_public_aux(name, child)
            end
        end

        return nothing
    end

    origin = eval(module_name)
    @assert origin isa Module

    for name in names(origin; all = true)
        make_public_aux(name, origin)
    end

    return nothing
end
export @make_public_rec

## function for threaded version of pmap (taken from https://github.com/JuliaLang/julia/issues/17887#issuecomment-564844185)
function tmap(f, xs::AbstractArray)
    g = Base.Generator(f,xs)
    et = Base.@default_eltype(g)
    a = Array{et}(undef, length(xs))
    Threads.@threads for i in 1:length(xs)
        a[i] = f(xs[i])
    end
    a
end

function permute_all(F::AbstractVector,m::Int;replace::Bool=false)
    ##find all sets of objects F choose m
    ## F is the set of objects, m is how many we choose for each permutation
    #if replace is false, m!>|F|
    if((replace == false) & (m>length(F)))
        throw(ArgumentError("When replace=false, m!>|F| (i.e. cannot choose more elements than provided)"))
    end
    n = length(F)
    comb = []
    for i in 1:m
        push!(comb,repeat(vcat([repeat([F[x]],n^(m-i)) for x in 1:n]...),n^(i  -1)))
    end
    candidates = hcat(comb...)
    if replace==false
        candidates = candidates[length.(unique.(eachrow(candidates))).==m,:]   
    end
    return candidates
end;   
export permute_all


function data_from_dataframe(df::DataFrame,identifier::String="data")
    #transform the numerical data columns of a DataFrame into an array. The identifier is a string that is common and uniquely contained in the names of the desired columns 
    data = Array(select(df,filter(x->occursin(identifier,x),names(df))))
end
export data_from_dataframe

end # module



