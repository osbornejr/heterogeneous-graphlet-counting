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
    if(src!="Initialisation.jl")
        includet("$cwd/src/"*src)
    end
end

#set up distributed workers
#first clean to make sure there are no stray workers already around
if(length(workers())!=Threads.nthreads())
    rmprocs(workers())
    #add workers equal to the number of available cpus  
    addprocs(Threads.nthreads())
    #addprocs(8)
    @everywhere include("src/CoexpressionMeasures.jl")
    @everywhere include("src/GraphletCounting.jl")
    @everywhere include("src/NullModel.jl")

end
