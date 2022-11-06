
types = ["coding","noncoding"]
colours = ["coding","noncoding"]
orbits = Dict{String,BitMatrix}("2-path"=>[1,1],
                                 "3-path"=>[1,2,1],
                                 "3-tri"=>[1,1,1],
                                 "4-path"=>[1,2,2,1],
                                 "4-star"=>[1,1,2,1],
                                 "4-tail"=>[1,1,2,3],
                                 "4-cycle"=>[1,2,1,2],
                                 "4-chord"=>[1,2,1,2],
                                 "4-clique"=>[1,1,1,1]
                                )
###note that the above orbit input for 4-cycle is technically incorrect, but it allows us to generate in this case the sole(?) het graphlet that is not implied by orbit structure. It also produces 3 other non-unique-symmetry 4-cycle graphlets that can be ignored.  
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
    ## colours need to be reversed if graphlet is all noncoding TODO currently only works if number of types is 2; make work for any number of types
    if ((!occursin("_coding",g)) & (!startswith(g,"coding")))
        NetworkConstruction.draw_tex_graphlet(g,colours=reverse(colours),out_file = out*g*".tex")
    else
        NetworkConstruction.draw_tex_graphlet(g,colours=colours,out_file = out*g*".tex")
    end
end

#save list to file as well
outfile= out*"list.txt"
open(outfile, "w") do f
  for i in List
    println(f, i)
  end
end 
