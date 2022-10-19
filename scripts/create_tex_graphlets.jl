
types = ["coding","noncoding"]
colours = ["coding","noncoding"]
orbits = Dict{String,Vector{Int}}("2-path"=>[1,1],
                                 "3-path"=>[1,2,1],
                                 "3-tri"=>[1,1,1],
                                 "4-path"=>[1,2,2,1],
                                 "4-star"=>[1,1,2,1],
                                 "4-tail"=>[1,1,2,3],
                                 "4-cycle"=>[1,2,1,2],
                                 "4-chord"=>[1,2,2,1],
                                 "4-clique"=>[1,1,1,1]
                                )
#generate list of all required graphlets
List = String[]
for (g,o) in orbits
    #generate list and append graphlet name
    list = NetworkConstruction.generate_heterogeneous_graphlet_list(o,types)
    list = list.*"_$(g)"
    for x in list
        push!(List,x)
    end
end
#where to store tex pics
out = ENV["PWD"]*"/output/share/graphlets/" 
run(`mkdir -p $(out)`)
for g in List
    NetworkConstruction.draw_tex_graphlet(g,colours=colours,out_file = out*g*".tex")
end
#save list to file as well
outfile= out*"list.txt"
open(outfile, "w") do f
  for i in List
    println(f, i)
  end
end 
