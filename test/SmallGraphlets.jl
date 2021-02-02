using Test
include("src/GraphletCounting.jl")

##A minimal test for each graphlet. Initially with two nodes of each colour
testvlist = [ "1" "red"; "2" "red";"3" "blue"; "4" "blue"]

## 4-tail test{{{
#     1-------2
#     | \
#     |   \
#     |     \
#     3-------4

testelist = [ 1=>2;1 =>3 ;1=>4;3=>4]

four_tail_expected = Dict{String,Int}("blue_blue_red_3-tri"=>1,"blue_blue_red_red_4-tail"=>1,"blue_red_red_3-path" =>2)

four_tail_test = count_graphlets(testvlist[:,2],testelist,4)
@test four_tail_test == four_tail_expected
#=}}}=#

## 4-star test{{{
#     1-------2
#     | \
#     |   \
#     |     \
#     3       4
	
testelist = [ 1=>2;1 =>3 ;1=>4]

four_star_expected = Dict{String,Int}("blue_blue_red_3-tri"=>1,"blue_blue_red_red_4-tail"=>1,"blue_red_red_3-path" =>2)

four_tail_test = count_graphlets(testvlist[:,2],testelist,4)
@test four_tail_test == four_tail_expected
#=}}}=#
