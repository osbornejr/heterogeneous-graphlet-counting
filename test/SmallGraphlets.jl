using Test
include("src/GraphletCounting.jl")

#### A minimal test for each graphlet. Initially just with two types (colours).

## 4-path tests{{{


#     1R-----2R
#     |  
#     |    
#     |      
#     3R-----4R

testvlist = [ "1" "red"; "2" "red";"3" "red"; "4" "red"]
testelist = [ 1=>2;1 =>3 ;3=>4]

four_path_expected_1 = Dict{String,Int}("red_red_red_red_4-path"=>1,"red_red_red_3-path" =>2)

four_path_test_1 = count_graphlets(testvlist[:,2],testelist,4)
#-----------------------------------------------------------------------------------------
#     1R-----2B
#     |  
#     |    
#     |      
#     3R-----4R

testvlist = [ "1" "red"; "2" "blue";"3" "red"; "4" "red"]
testelist = [ 1=>2;1 =>3 ;3=>4]

four_path_expected_2 = Dict{String,Int}("blue_red_red_red_4-path"=>1,"red_red_red_3-path" =>1,"blue_red_red_3-path" =>1)

four_path_test_2 = count_graphlets(testvlist[:,2],testelist,4)
#--------------------------------------------------------------------------------------------

#     1B-----2R
#     |  
#     |    
#     |      
#     3R-----4R

testvlist = [ "1" "blue"; "2" "red";"3" "red"; "4" "red"]
testelist = [ 1=>2;1 =>3 ;3=>4]

four_path_expected_3 = Dict{String,Int}("red_blue_red_red_4-path"=>1,"red_blue_red_3-path" =>1,"blue_red_red_3-path" =>1)

four_path_test_3 = count_graphlets(testvlist[:,2],testelist,4)
#--------------------------------------------------------------------------------------------

#     1B-----2B
#     |  
#     |    
#     |      
#     3R-----4R

testvlist = [ "1" "blue"; "2" "blue";"3" "red"; "4" "red"]
testelist = [ 1=>2;1 =>3 ;3=>4]

four_path_expected_4 = Dict{String,Int}("blue_blue_red_red_4-path"=>1,"blue_blue_red_3-path" =>1,"blue_red_red_3-path" =>1)

four_path_test_4 = count_graphlets(testvlist[:,2],testelist,4)

#--------------------------------------------------------------------------------------------

#     1R-----2B
#     |  
#     |    
#     |      
#     3B-----4R

testvlist = [ "1" "red"; "2" "blue";"3" "blue"; "4" "red"]
testelist = [ 1=>2;1 =>3 ;3=>4]

four_path_expected_5 = Dict{String,Int}("red_blue_red_blue_4-path"=>1,"red_blue_red_3-path" =>1,"blue_red_blue_3-path" =>1)

four_path_test_5 = count_graphlets(testvlist[:,2],testelist,4)

#--------------------------------------------------------------------------------------------
#     1R-----2B
#     |  
#     |    
#     |      
#     3R-----4B

testvlist = [ "1" "red"; "2" "blue";"3" "red"; "4" "blue"]
testelist = [ 1=>2;1 =>3 ;3=>4]

four_path_expected_6 = Dict{String,Int}("blue_red_red_blue_4-path"=>1,"blue_red_red_3-path" =>2)

four_path_test_6 = count_graphlets(testvlist[:,2],testelist,4)
#--------------------------------------------------------------------------------------------

@test four_path_test_1 == four_path_expected_1
@test four_path_test_2 == four_path_expected_2
@test four_path_test_3 == four_path_expected_3
@test four_path_test_4 == four_path_expected_4
@test four_path_test_5 == four_path_expected_5
@test four_path_test_6 == four_path_expected_6#=}}}=#

## 4-tail tests{{{

#     1R-----2R
#     | \
#     |   \
#     |     \
#     3R-----4R

testvlist = [ "1" "red"; "2" "red";"3" "red"; "4" "red"]
testelist = [ 1=>2;1 =>3 ;1=>4;3=>4]

four_tail_expected_1 = Dict{String,Int}("red_red_red_3-tri"=>1,"red_red_red_red_4-tail"=>1,"red_red_red_3-path" =>2)

four_tail_test_1 = count_graphlets(testvlist[:,2],testelist,4)


#     1R-----2B
#     | \
#     |   \
#     |     \
#     3R-----4R

testvlist = [ "1" "red"; "2" "blue";"3" "red"; "4" "red"]
testelist = [ 1=>2;1 =>3 ;1=>4;3=>4]

four_tail_expected_2 = Dict{String,Int}("red_red_red_3-tri"=>1,"red_red_red_blue_4-tail"=>1,"blue_red_red_3-path" =>2)

four_tail_test_2 = count_graphlets(testvlist[:,2],testelist,4)

#     1B-----2R
#     | \
#     |   \
#     |     \
#     3R-----4R

testvlist = [ "1" "blue"; "2" "red";"3" "red"; "4" "red"]
testelist = [ 1=>2;1 =>3 ;1=>4;3=>4]

four_tail_expected_3 = Dict{String,Int}("blue_red_red_3-tri"=>1,"red_red_blue_red_4-tail"=>1,"red_blue_red_3-path" =>2)

four_tail_test_3 = count_graphlets(testvlist[:,2],testelist,4)

#     1R-----2R
#     | \
#     |   \
#     |     \
#     3B-----4R

testvlist = [ "1" "red"; "2" "red";"3" "blue"; "4" "red"]
testelist = [ 1=>2;1 =>3 ;1=>4;3=>4]

four_tail_expected_4 = Dict{String,Int}("blue_red_red_3-tri"=>1,"blue_red_red_red_4-tail"=>1,"blue_red_red_3-path" =>1,"red_red_red_3-path"=>1)

four_tail_test_4 = count_graphlets(testvlist[:,2],testelist,4)

#     1B-----2B
#     | \
#     |   \
#     |     \
#     3R-----4R

testvlist = [ "1" "blue"; "2" "blue";"3" "red"; "4" "red"]
testelist = [ 1=>2;1 =>3 ;1=>4;3=>4]

four_tail_expected_5 = Dict{String,Int}("blue_red_red_3-tri"=>1,"red_red_blue_blue_4-tail"=>1,"blue_blue_red_3-path" =>2)

four_tail_test_5 = count_graphlets(testvlist[:,2],testelist,4)

#     1R-----2B
#     | \
#     |   \
#     |     \
#     3B-----4R

testvlist = [ "1" "red"; "2" "blue";"3" "blue"; "4" "red"]
testelist = [ 1=>2;1 =>3 ;1=>4;3=>4]

four_tail_expected_6 = Dict{String,Int}("blue_red_red_3-tri"=>1,"blue_red_red_blue_4-tail"=>1,"blue_red_blue_3-path" =>1,"blue_red_red_3-path"=>1)

four_tail_test_6 = count_graphlets(testvlist[:,2],testelist,4)

#     1B-----2R
#     | \
#     |   \
#     |     \
#     3B-----4R

testvlist = [ "1" "blue"; "2" "red";"3" "blue"; "4" "red"]
testelist = [ 1=>2;1 =>3 ;1=>4;3=>4]

four_tail_expected_7 = Dict{String,Int}("blue_blue_red_3-tri"=>1,"blue_red_blue_red_4-tail"=>1,"blue_blue_red_3-path" =>1,"red_blue_red_3-path"=>1)

four_tail_test_7 = count_graphlets(testvlist[:,2],testelist,4)

#     1R-----2R
#     | \
#     |   \
#     |     \
#     3B-----4B

testvlist = [ "1" "red"; "2" "red";"3" "blue"; "4" "blue"]
testelist = [ 1=>2;1 =>3 ;1=>4;3=>4]

four_tail_expected_8 = Dict{String,Int}("blue_blue_red_3-tri"=>1,"blue_blue_red_red_4-tail"=>1,"blue_red_red_3-path" =>2)

four_tail_test_8 = count_graphlets(testvlist[:,2],testelist,4)

@test four_tail_test_1 == four_tail_expected_1
@test four_tail_test_2 == four_tail_expected_2
@test four_tail_test_3 == four_tail_expected_3
@test four_tail_test_4 == four_tail_expected_4
@test four_tail_test_5 == four_tail_expected_5
@test four_tail_test_6 == four_tail_expected_6
@test four_tail_test_7 == four_tail_expected_7
@test four_tail_test_8 == four_tail_expected_8
#=}}}=#

## 4-star test{{{
#     1R-----2R
#     | \
#     |   \
#     |     \
#     3B     4B
    
testelist = [ 1=>2;1 =>3 ;1=>4]

four_star_expected = Dict{String,Int}("blue_red_blue_3-path"=>1,"blue_blue_red_red_4-star"=>1,"blue_red_red_3-path" =>2)

four_star_test = count_graphlets(testvlist[:,2],testelist,4)
@test four_star_test == four_star_expected
#=}}}=#

#Test a full run of the old method against the new method   
g = Vector{String}()#={{{=#
d =  Vector{Int}()
for graphlet in collect(keys(graphlet_counts))
    if graphlet_counts[graphlet]!==graphlet_counts_old[graphlet]
        push!(g,graphlet)
        push!(d,graphlet_counts[graphlet]-graphlet_counts_old[graphlet])
    end
end
errors = hcat(g,d)
@test length(errors)==0#=}}}=#

#Test a full run of the single thread method against multi thread method    
g = Vector{String}()#={{{=#
d =  Vector{Int}()
for graphlet in collect(keys(graphlet_counts_single))
    if graphlet_counts_single[graphlet]!==graphlet_counts_threaded[graphlet]
        push!(g,graphlet)
        push!(d,graphlet_counts_single[graphlet]-graphlet_counts_threaded[graphlet])
    end
end
errors = hcat(g,d)
@test length(errors)==0#=}}}=#

