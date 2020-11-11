
#### RT-POL specific archive
##import data
import CSV
edgelist=CSV.read("data/rt-pol/soc-political-retweet.edges",header=["V1","V2","ref"],type=Int)
edgelist=edgelist[:,1:2]
vertexlist=CSV.read("data/rt-pol/soc-political-retweet.node_labels",header=["V","type"],type=Int)

#fix vertex values to iterate from 1, not 0
edgelist[:,1:2]=edgelist[:,1:2].+1
vertexlist.V=vertexlist.V.+1

#change label names
vertexlist.type=replace(vertexlist.type,Pair(1,"right"),Pair(2,"left"))


#start writing required functions ad hoc
#change to arrays for performance (also remove multi-edges from edgelist here for now)
edgelist=unique(sort(convert(Array{Int,2},edgelist),dims=2),dims=1)
#we now need the edgelist to consist of pairs
splat(f) = args->f(args...)
edgelist = splat(Pair).(eachrow(edgelist))


orbit_counts = count_graphlets(String.(vertexlist[:,2]),edgelist,3)
#count triangles for comparison
tri_counts=[orbit_counts[i] for i in collect(keys(orbit_counts)) if occursin("3-tri",i)]
tri_counts=tri_counts./sum(tri_counts)
