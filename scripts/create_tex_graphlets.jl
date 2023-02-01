
types = ["coding","noncoding"]
colours = ["coding","noncoding"]

#Input the adj matrix (upper triangle fill style) for each graphlet. Note, this must match the configuration baked into the draw_tex_graphlets function of
#          1--------4
#          |        |        
#          |        |        
#          |        |        
#          2--------3
adjs = Dict{String,BitVector}("2-path"=>[true],
                                 "3-path"=>[true,false,true],
                                 "3-tri"=>[true,true,true],
                                 "4-path"=>[true,false,false,true,false,true],
                                 "4-star"=>[false,true,false,true,false,true],
                                 "4-tail"=>[true,true,false,true,false,true],
                                 "4-cycle"=>[true,false,true,true,false,true],
                                 "4-chord"=>[true,true,true,true,false,true],
                                 "4-clique"=>[true,true,true,true,true,true]
                                )
#generate list of all required graphlets
List = String[]
for (g,o) in adjs
    #generate list and append graphlet name
    list = NetworkConstruction.generate_heterogeneous_graphlet_list(o,types)
    list = list.*"_$(g)"
    for x in list
        push!(List,x)
    end
end


#where to store tex pics
out = ENV["PWD"]*"/output/share/graphlets/" 
#clear and remake directory to update
run(`rm -r $(out)`)
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
