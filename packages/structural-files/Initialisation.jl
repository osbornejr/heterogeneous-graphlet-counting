#need to activate Revise here for now (because Pluto?)
using Revise
### Include all source files TODO make this occur more fluently and automatically by creating a package, and using Revise
cwd = ENV["JULIA_PROJECT"]


for src in filter(x->endswith(x,".jl"),readdir("src"))
    if(!occursin("Initialisation.jl",src))
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
