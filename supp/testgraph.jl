using Plots, GraphPlot,  LightGraphs, Compose, Measures, Printf
import Cairo
anim=Animation()
g = graphfamous("karate")
l_x = 2 * rand(nv(g)) .- 1.0
l_y = 2 * rand(nv(g)) .- 1.0
for i in 1:10
    mylayout(g, kws...) = spring_layout(g, l_x, l_y, kws...)
    p=gplot(g, layout=mylayout, 
        nodelabel = 1:nv(g))
    output = compose(p,
        (context(), Compose.text(1, 1, "Julia")),
        (context(), rectangle(), fill("white")))
    j=length(anim.frames) + 1
    tmpfilename=joinpath(anim.dir,@sprintf("%06d.png",j))
    Compose.draw(PNG(tmpfilename),output)
    push!(anim.frames, tmpfilename)
end
gif(anim, "anim_test.gif", fps = 5)
