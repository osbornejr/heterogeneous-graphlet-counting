### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ d4f1c85c-f854-11ed-09a1-af06379be62d
# ╠═╡ show_logs = false
begin
cwd = ENV["PWD"];	
	import Pkg
	#Pkg.activate(mktempdir())
    #Pkg.add([ 
    #    Pkg.PackageSpec(name="DataFrames",version'"1.3"),
	#	Pkg.PackageSpec(name="Revise"),
	#	Pkg.PackageSpec(name="CommonMark"),
	#	Pkg.PackageSpec(name="PlutoUI"),
	#Pkg.PackageSpec(name="DataStructures"),
	#	Pkg.PackageSpec(name="StatsBase"),
	#	Pkg.PackageSpec(name="Gadfly"),
	#	Pkg.PackageSpec(name="LinearAlgebra")
    #])
	using Revise	
	using CommonMark
	using PlutoUI
	using LinearAlgebra
	using DataFrames
	using DataStructures
	using StatsBase
	using Gadfly
	using Portinari
	using HypertextLiteral: JavaScript,@htl,@htl_str
	## load dev packages last
	dev = Pkg.develop
	dev(path=cwd*"/packages/DataPreprocessing")
	dev(path=cwd*"/packages/NetworkConstruction")
	dev(path=cwd*"/packages/GraphletCounting")
	dev(path=cwd*"/packages/GraphletAnalysis")
	dev(path=cwd*"/packages/ProjectFunctions")
	using DataPreprocessing
	using NetworkConstruction
	using GraphletCounting
	using GraphletAnalysis
	using ProjectFunctions
	
end

# ╔═╡ 9752572a-cc9b-421d-ba41-05a1a8b6def6
md"""
## Parameters
"""

# ╔═╡ ee05e4b5-13f4-433a-9076-ba304c5c8807
md"""
### Experiment
"""

# ╔═╡ 32668349-8283-40c6-bf0e-1afc5ee960ea
@bind experiment Select(["GSE68559","GSE68559_sub"])

# ╔═╡ 351ea0cc-ce6b-4e82-99dd-8773c5a9e8a4
load_config(cwd*"/config/run-files/$experiment.yaml")

# ╔═╡ aecb1084-001e-4b0d-a3f0-740de58e5753
md"""
### Threshold on per sample significance
"""

# ╔═╡ 07904d8d-2e38-4207-aa94-0dd835708854
@bind expression_cut Select([0.05,0.10])

# ╔═╡ 52ca87b1-8796-4c0c-9d8c-95532d226e52
md"""
### Minimum samples required
"""

# ╔═╡ d9a4d419-ed7e-4dc9-aa50-d91c166a77d8
@bind minreq Select([0.05,0.10])

# ╔═╡ b6c9b8a2-3b94-4570-adff-c3ae0962683d
md"""
### Normalisation method
"""

# ╔═╡ 98c53468-1da4-4515-a135-5445eba725b5
@bind norm_meth Select(["median","upper_quartile","quantile","TMM","TMMwsp","total_count"])

# ╔═╡ 4ea93937-05ac-4f2e-bd66-b6898489c24d
md"""
### Variance cut off
"""

# ╔═╡ 57477b4e-2a67-4817-8ef3-247c6ba2f8e9
@bind variance_cut Select([1.0,2.0])

# ╔═╡ 4aecc101-7800-40c0-a8b0-8e0f67e98cb3
# ╠═╡ show_logs = false
begin
	params["data_preprocessing"]["expression_cutoff"] = expression_cut
	params["data_preprocessing"]["minreq"] = minreq
	params["data_preprocessing"]["norm_method"] = norm_meth
	params["data_preprocessing"]["variance_percent"] = variance_cut 
	ProjectFunctions.cache_setup()
	raw_counts,round_counts,vst_counts,clean_counts,norm_counts,processed_counts = get_preprocessed_data();
end;

# ╔═╡ 20215a75-b329-48fa-9024-bd6aab7f8247
params

# ╔═╡ ee390e38-f48f-4ab9-9947-87de4224ca94
summarystats(raw_counts."GSM1675513_MB011_1 data")

# ╔═╡ 84bc738a-4e89-4811-88b6-ec12b81bed1f
test = processed_counts

# ╔═╡ 0686222b-0955-4421-b727-e9c57b9a352d
sum(raw_counts.transcript_type.=="noncoding")

# ╔═╡ 104f6b63-70e0-424f-b549-84c4d8334acb
cod = sum(test.transcript_type.=="coding")

# ╔═╡ b1984144-dbbc-4a00-9273-c4ffa76d2f53
nod = sum(test.transcript_type.=="noncoding")

# ╔═╡ 3b942054-f8ca-4449-a21d-a0d2e984b0ad
nod/(nod+cod)

# ╔═╡ 8c563073-fb2d-40e5-9113-829f63594205
DataPreprocessing.boxplot(raw_counts);

# ╔═╡ 739173a4-5f67-416b-b276-774c21aae4f9
DataPreprocessing.boxplot(clean_counts);

# ╔═╡ f8c20e2e-252c-4796-8dd6-c910eb6f7820
DataPreprocessing.boxplot(norm_counts);

# ╔═╡ b4373d9d-e7da-4fee-a487-1c5ceae54c61
    regions= [
		  filter(x->occursin("data",x),names(raw_counts))=>"All",
      filter(x->occursin("_1 ",x),names(raw_counts))=>"BA10",
	filter(x->occursin("_2 ",x),names(raw_counts))=>"BA22",
	filter(x->occursin("_3 ",x),names(raw_counts))=>"BA24",
	filter(x->occursin("_4 ",x),names(raw_counts))=>"insula",
	filter(x->occursin("_5 ",x),names(raw_counts))=>"amygdala",
	filter(x->occursin("_6 ",x),names(raw_counts))=>"hippocampus",
	filter(x->occursin("_7 ",x),names(raw_counts))=>"posterior putamen",
	filter(x->occursin("_8 ",x),names(raw_counts))=>"cerebellum",
	filter(x->occursin("_9 ",x),names(raw_counts))=>"raphae nuclei",
	filter(x->occursin("_10 ",x),names(raw_counts))=>"BA46"];


# ╔═╡ c7fadafb-45c5-4168-8310-4d6dd5fcdd38
region_dict = Dict(Pair.(last.(regions),first.(regions)));
	

# ╔═╡ 2b050ca8-3050-4e8a-bd00-8ab9b40001ed
    patients = [
		filter(x->occursin("data",x),names(raw_counts))=>"All",
	 filter(x->occursin("MB011",x),names(raw_counts))=>"smoker MB011" ,
	 filter(x->occursin("MB059",x),names(raw_counts))=>"smoker MB059" ,
	 filter(x->occursin("MB100",x),names(raw_counts))=>"smoker MB100" ,
	 filter(x->occursin("MB148",x),names(raw_counts))=>"smoker MB148" ,
	 filter(x->occursin("MB160",x),names(raw_counts))=>"smoker MB160" ,
	 filter(x->occursin("MB052",x),names(raw_counts))=>"nonsmoker MB052" ,
	 filter(x->occursin("MB147",x),names(raw_counts))=>"nonsmoker MB147" ,
	 filter(x->occursin("MB151",x),names(raw_counts))=>"nonsmoker MB151" ,
	 filter(x->occursin("MB197",x),names(raw_counts))=>"nonsmoker MB197" ,
	 filter(x->occursin("MB202",x),names(raw_counts))=>"nonsmoker MB202"];
	

# ╔═╡ 6dcf8d2a-1f96-46a8-bf2c-57feea99dad3
patient_dict = Dict(Pair.(last.(patients),first.(patients)));

# ╔═╡ 3e05e27d-0445-4773-b242-bb3e7a4853ca
@bind REG Select(regions)

# ╔═╡ ded8aa8b-db1a-4f2a-abe0-b3a8d3cc7ca1
@bind PAT Select(patients)

# ╔═╡ 2beaf388-3d97-41da-ba16-41154906155d
selected_names = intersect(REG,PAT)

# ╔═╡ 5e7acfb3-4566-4a65-9798-85bf0584b439
#titles = ["VST counts","Cleaned counts", "Normalised counts ($norm_meth)"]
titles = ["Raw counts","Rounded counts","VST counts","Cleaned counts", "Normalised counts ($norm_meth)","Sampled counts"]


# ╔═╡ d9c0a312-4b19-49c7-8f0f-2d1fcfa13b31
@htl("""
		<style>	.plutoui-rangeslider { width: 50em } </style>
		<h3 class="dash">$(titles[index])</h3>
		<div style="display: flex; justify-content: center; align-items: center; gap: 2em"></div>
	""")

# ╔═╡ b1fba899-0ca7-48ba-bbe0-f30a146536a1
@htl(
	"""
<meta charset="utf-8">

<!-- Load d3.js -->
<script src="https://d3js.org/d3.v4.js"></script>
 
<!-- Create a div where the graph will take place -->
<div id="my_dataviz"></div>
<script> 
$(JavaScript(text)) 
</script>	
""" 
)

# ╔═╡ fbdd59ff-9b6c-4498-bdb7-b36682b45ed5
@bind tick Clock()

# ╔═╡ 8e7533fa-4914-44f9-b470-fea060af8fe1
@bind step Select([1=>"raw counts",2=>"round counts",3=>"vst counts",4=>"clean counts",5=>"norm counts",6=>"processed counts"])

# ╔═╡ d4ab62e9-de6e-42fa-a655-d2c8231b15c6
@bind index_switch Select([0=>"select",1=>"clock"])

# ╔═╡ f2eb27c9-e0e1-4687-aed8-927df4d5decc
begin
	#set = [vst_counts,clean_counts,norm_counts]
	set = [raw_counts,round_counts,vst_counts,clean_counts,norm_counts,processed_counts]
	index = index_switch*(mod(tick,length(set))+1) + (-1*index_switch+1)*(step)
	input = set[index]
	width = min(1600,400*length(selected_names))
	height = 400
	range = max(data_from_dataframe(input)...)	
	text = "// set the dimensions and margins of the graph 
var margin = {top: 10, right: 30, bottom: 30, left: 40}, 
  width = $(width) - margin.left - margin.right, 
  height =$(height) - margin.top - margin.bottom; 
 
// append the svg object to the body of the page 
var svg = d3.select(\"#my_dataviz\") 
.append(\"svg\") 
  .attr(\"width\", width + margin.left + margin.right) 
  .attr(\"height\", height + margin.top + margin.bottom) 
.append(\"g\") 
  .attr(\"transform\", 
        
	\"translate(\" + margin.left + \",\" + margin.top + \")\"); 
 "
	
	for (i,name) in enumerate(selected_names)
		#i = 1
		#name = names(raw_counts)[5]
		#sample data (in FPKM)
		data = input[:,name]
		stats = summarystats(data)
		
		text*="
// create dummy data 
// var data = [12,3,4,4] 
  
// Compute summary statistics used for the box: 
var q1 = $(stats.q25) 
var median = $(stats.median) 
var q3 = $(stats.q75) 
var interQuantileRange = q3 - q1 
var min = $(max(stats.min,stats.q25-1.5*(stats.q75-stats.q25))) 
var max = $(min(stats.max,stats.q75+1.5*(stats.q75-stats.q25)))  
 
// Show the Y scale 
var y = d3.scaleLinear() 
		.domain([$(stats.min),4]) 
  .range([height, 0]); 
var axis = svg.call(d3.axisLeft(y)) 
 
  
  axis.selectAll(\"line\") 
    .style(\"stroke\", \"grey\"); 
 
  axis.selectAll(\"path\") 
    .style(\"stroke\", \"grey\"); 
 
  axis.selectAll(\"text\") 
    .style(\"stroke\", \"grey\"); 
	 
// a few features for the box 
var center = $(i*(width-70)/length(selected_names)) 
var width = $((width-70)/length(selected_names)-0.2*(width-70)/length(selected_names)) 
 
// Show the main vertical line 
svg 
.append(\"line\") 
  .attr(\"x1\", center) 
  .attr(\"x2\", center) 
  .attr(\"y1\", y(min) ) 
  .attr(\"y2\", y(max) ) 
  .attr(\"stroke\", \"grey\") 
 
// Show the box 
svg 
.append(\"rect\") 
  .attr(\"x\", center - width/2) 
  .attr(\"y\", y(q3) ) 
  .attr(\"height\", (y(q1)-y(q3)) ) 
  .attr(\"width\", width ) 
.style(\"fill\", \"#69e3a2\") 

 
 
// show median, min and max horizontal lines 
svg 
.selectAll(\"toto\") 
.data([min, median, max]) 
.enter() 
.append(\"line\") 
  .attr(\"x1\", center-width/2) 
  .attr(\"x2\", center+width/2) 
  .attr(\"y1\", function(d){ return(y(d))} ) 
  .attr(\"y2\", function(d){ return(y(d))} ) 
  .attr(\"stroke\", \"grey\") 
" 
	end
end

# ╔═╡ 332cf792-e053-41d3-a1f0-cb0ff90311b7
md"""
## $(titles[index])
"""

# ╔═╡ b4262a92-7bc8-44a3-acaf-3dbcc235f688
sum(data_from_dataframe(raw_counts),dims=1)

# ╔═╡ 4b63bfcf-d192-4826-ab41-84b0cc1b6489
raw_counts

# ╔═╡ 757bb449-28ec-4545-b1b5-d3f84217ab81
begin
	index
	TableOfContents()
end

# ╔═╡ 92c7a3a4-ebc1-4166-a106-9b781f546f60
vec(sum(DataPreprocessing.data_from_dataframe(norm_counts).==0,dims=1))

# ╔═╡ ffe2cd5b-4be6-40c7-aa97-72db54fd9d3a
sum(input.transcript_type.=="coding")/size(input)[1]

# ╔═╡ 1ca33a55-85ab-47ef-9966-64fcc963c0be
input_data = data_from_dataframe(input,"data");

# ╔═╡ 3557b169-6656-49a5-b9eb-4f3d0ced356a
sum(input_data.==0.0)/length(input_data)

# ╔═╡ 048bb985-d7b8-4ace-8606-8ddfddc51ba1
begin
	#test = DataPreprocessing.clean_raw_counts(new_clean_counts,1)
	data = filter(>(0),input[:,2])
	h = fit(Histogram,data,nbins=100)
end

# ╔═╡ cfe7529a-ca49-4329-846d-9e4cfcd43e8b
h.edges

# ╔═╡ 7a2168ad-61a0-4132-90ea-69574f0040fc
Bars(round.(collect(first(h.edges))[2:end],digits = 3).|>string,h.weights,"histogram";attributes=D3Attr(attr=(;fill="rgba(10, 200, 100, 0.6)")))

# ╔═╡ 9863d522-bb4f-4287-8c7b-06461be85369
sig = params["data_preprocessing"]["roundsig"]

# ╔═╡ 46085b1f-a043-45ad-96db-ae1a8e7dceea
md"""
## Rounding raw data
To increase clarity and avoid skewing towards variance of lowly expressed transcripts, we round the raw FPKM values to $(sig) significant figures.
The main impact of this rounding at this stage is to define which expression entries are considered nonzero i.e. any expression reading ``\geq`` $(10.0^(-1*sig)) is considered nonzero. 
"""




# ╔═╡ 3d55ece3-b8da-4443-91c9-98b9b983a04e
md"""
## Variance stabilising transformation
We apply a variance stabilising transformation. 
In this case this is a simple log-transform $$\log_2(x+1)$$.
"""

# ╔═╡ 692edabc-3542-49fa-b532-c5172f039951
md"""
## Clean expression data
Our next task is to remove any transcripts that have a low expression profile across samples.
To do this, we select a cut off point that represents the expression value that is the lower bound for the top 10% expression values across ALL non-zero entries in all samples.
We then compare each transcript's expression profile to this cutoff value.
If at least 10% of its expression values are greater than the cutoff value, we keep the transcript, but otherwise we discard.
"""

# ╔═╡ e2c4ec1b-bd5f-4d49-a6cc-b9366cfdc9c5
md"""
## Normalising data

One question here is whether it matters if we are normalising VST or non VST data here.

Will the difference matter?

Following the basic model used in `edgeR` normalisation, we normalise using library size 
"""

# ╔═╡ c8cedcba-6d09-40c7-905c-a6a47596798d
md"""
## PCA check

"""

# ╔═╡ 704d588e-cffd-4239-b5c9-2760d8afdb6f
PCs,D,per_sample=DataPreprocessing.pca(data_from_dataframe(input));

# ╔═╡ 87b0e404-4118-4abc-9f2e-4a963bd4496d
sample_names = filter(x->occursin("data",x),names(raw_counts))


# ╔═╡ 5bd9c340-4d63-4a54-be4e-8c121bebabd5
marg = 10

# ╔═╡ 818cadda-b311-4aae-9c33-2d85a18899f5
a = 1

# ╔═╡ 97235b4e-93fc-4ed4-a3a2-c6541cb9ef68
b = 2

# ╔═╡ 4e3eadbe-7960-42e9-b399-3a44297d70f0
nx = collect(min(per_sample[:,b]...)-marg:max(per_sample[:,b]...)+marg);

# ╔═╡ 1c0643d0-332f-45d5-b7e8-0b9b6b6b680b
ny = collect(min(per_sample[:,a]...)-marg:max(per_sample[:,a]...)+marg)

# ╔═╡ 2ffc3d13-ab1b-4b03-864a-66a0896dbdff
colour = ["114,245,51","64,165,245","245,223,38","245,93,239","245,159,27","245,154,51","245,250,64","245,40,82","71,195,245","101,27,245"]

# ╔═╡ c211987a-bc3f-4ade-a5a0-e19563cf9ddd
findall(x->x==min(per_sample[:,3]...),per_sample[:,3])

# ╔═╡ 75641449-1801-4f22-a90b-5947f6c34b8a
sample_names[7]

# ╔═╡ 967ffab4-8abf-4447-8f09-cd0e5d624767
norm_counts

# ╔═╡ 01ff1c3e-6eaa-4863-bb1d-7f6f7037fafe
begin
	var_data = vec(var(data_from_dataframe(norm_counts),dims =2 ))
	v = fit(Histogram,var_data,nbins=100)
	Bars(round.(collect(first(v.edges))[2:end],digits = 3).|>string,v.weights,"histogram";attributes=D3Attr(attr=(;fill="rgba(10, 200, 100, 0.6)")))
end

# ╔═╡ d3762b94-02d7-466c-be63-a534d6a75894
processed_counts

# ╔═╡ 1414f977-cb37-4592-a104-217822a42ae2
v.edges

# ╔═╡ c0213d48-97da-47e3-a32c-9bf87c465c09
1600*1599/2

# ╔═╡ 07f12515-a320-434c-b1e3-73a305491f94
v.weights

# ╔═╡ e5e63932-1f89-4858-ae76-72fefa2d0139
norm_counts[var_data.>1.0,:]

# ╔═╡ 317602a3-d0b3-435f-af6d-112be42617d4
patient_vec =replace.(collect(keys(patient_dict))[2:end],"smoker "=>"","nonsmoker "=>"")

# ╔═╡ 21aa68b4-6257-41e6-a540-1f4bf92ef96c
Context(
	(;domain=([extrema(nx)...]), range=[0, 300]),
	(;domain=([extrema(ny)...]), range=[0, 300]),
[
Portinari.Shape(per_sample[findall(y->occursin("_$(patient_vec[i])",y),sample_names),b],per_sample[findall(y->occursin("_$(patient_vec[i])",y),sample_names),a],repeat([100],length(findall(y->occursin("_$(patient_vec[i])",y),sample_names))),"PC";
attributes=D3Attr(style=(;fill="rgba($(colour[i]), 0.8)"))
) for i in 1:10],
	"together"
)
	


# ╔═╡ a7ee49b3-3096-4818-8137-39ddc1d8195e
region_vec = collect(keys(region_dict))

# ╔═╡ 86e018ba-a78c-4644-af9e-63c6fec86e40
Portinari.Shape(per_sample[:,1],per_sample[:,3],repeat([100],length(sample_names)),"PC";
attributes=D3Attr(style=(;fill="rgba(255, 0, 0, 0.3)"))
)

# ╔═╡ a51fc8b8-e2d3-4962-abd9-54a5f07d517f
findall(x->x==max(abs.(per_sample[:,3])...),abs.(per_sample[:,3]))

# ╔═╡ c05fca81-7f29-4576-99d6-fc35210e25b7


# ╔═╡ d082290c-4412-42b5-aa15-e35e01923e5e
p_1=DataPreprocessing.pca_plot(per_sample,3)

# ╔═╡ 50958e94-9de4-4b76-8d79-3b7f65b94729
md"""
## Select for variation
"""

# ╔═╡ Cell order:
# ╠═d4f1c85c-f854-11ed-09a1-af06379be62d
# ╠═9752572a-cc9b-421d-ba41-05a1a8b6def6
# ╠═ee05e4b5-13f4-433a-9076-ba304c5c8807
# ╠═32668349-8283-40c6-bf0e-1afc5ee960ea
# ╠═351ea0cc-ce6b-4e82-99dd-8773c5a9e8a4
# ╠═aecb1084-001e-4b0d-a3f0-740de58e5753
# ╠═07904d8d-2e38-4207-aa94-0dd835708854
# ╠═52ca87b1-8796-4c0c-9d8c-95532d226e52
# ╠═d9a4d419-ed7e-4dc9-aa50-d91c166a77d8
# ╠═b6c9b8a2-3b94-4570-adff-c3ae0962683d
# ╠═98c53468-1da4-4515-a135-5445eba725b5
# ╠═4ea93937-05ac-4f2e-bd66-b6898489c24d
# ╠═57477b4e-2a67-4817-8ef3-247c6ba2f8e9
# ╠═4aecc101-7800-40c0-a8b0-8e0f67e98cb3
# ╠═20215a75-b329-48fa-9024-bd6aab7f8247
# ╠═ee390e38-f48f-4ab9-9947-87de4224ca94
# ╠═84bc738a-4e89-4811-88b6-ec12b81bed1f
# ╠═0686222b-0955-4421-b727-e9c57b9a352d
# ╠═104f6b63-70e0-424f-b549-84c4d8334acb
# ╠═b1984144-dbbc-4a00-9273-c4ffa76d2f53
# ╠═3b942054-f8ca-4449-a21d-a0d2e984b0ad
# ╠═8c563073-fb2d-40e5-9113-829f63594205
# ╠═739173a4-5f67-416b-b276-774c21aae4f9
# ╠═f8c20e2e-252c-4796-8dd6-c910eb6f7820
# ╠═c7fadafb-45c5-4168-8310-4d6dd5fcdd38
# ╠═b4373d9d-e7da-4fee-a487-1c5ceae54c61
# ╠═2b050ca8-3050-4e8a-bd00-8ab9b40001ed
# ╠═6dcf8d2a-1f96-46a8-bf2c-57feea99dad3
# ╠═3e05e27d-0445-4773-b242-bb3e7a4853ca
# ╠═ded8aa8b-db1a-4f2a-abe0-b3a8d3cc7ca1
# ╠═2beaf388-3d97-41da-ba16-41154906155d
# ╠═f2eb27c9-e0e1-4687-aed8-927df4d5decc
# ╠═5e7acfb3-4566-4a65-9798-85bf0584b439
# ╠═332cf792-e053-41d3-a1f0-cb0ff90311b7
# ╠═d9c0a312-4b19-49c7-8f0f-2d1fcfa13b31
# ╠═b1fba899-0ca7-48ba-bbe0-f30a146536a1
# ╠═fbdd59ff-9b6c-4498-bdb7-b36682b45ed5
# ╠═8e7533fa-4914-44f9-b470-fea060af8fe1
# ╠═d4ab62e9-de6e-42fa-a655-d2c8231b15c6
# ╠═b4262a92-7bc8-44a3-acaf-3dbcc235f688
# ╠═4b63bfcf-d192-4826-ab41-84b0cc1b6489
# ╠═757bb449-28ec-4545-b1b5-d3f84217ab81
# ╠═92c7a3a4-ebc1-4166-a106-9b781f546f60
# ╠═ffe2cd5b-4be6-40c7-aa97-72db54fd9d3a
# ╠═1ca33a55-85ab-47ef-9966-64fcc963c0be
# ╠═3557b169-6656-49a5-b9eb-4f3d0ced356a
# ╠═048bb985-d7b8-4ace-8606-8ddfddc51ba1
# ╠═cfe7529a-ca49-4329-846d-9e4cfcd43e8b
# ╠═7a2168ad-61a0-4132-90ea-69574f0040fc
# ╠═9863d522-bb4f-4287-8c7b-06461be85369
# ╠═46085b1f-a043-45ad-96db-ae1a8e7dceea
# ╠═3d55ece3-b8da-4443-91c9-98b9b983a04e
# ╠═692edabc-3542-49fa-b532-c5172f039951
# ╠═e2c4ec1b-bd5f-4d49-a6cc-b9366cfdc9c5
# ╠═c8cedcba-6d09-40c7-905c-a6a47596798d
# ╠═704d588e-cffd-4239-b5c9-2760d8afdb6f
# ╠═87b0e404-4118-4abc-9f2e-4a963bd4496d
# ╠═5bd9c340-4d63-4a54-be4e-8c121bebabd5
# ╠═4e3eadbe-7960-42e9-b399-3a44297d70f0
# ╠═818cadda-b311-4aae-9c33-2d85a18899f5
# ╠═97235b4e-93fc-4ed4-a3a2-c6541cb9ef68
# ╠═1c0643d0-332f-45d5-b7e8-0b9b6b6b680b
# ╠═21aa68b4-6257-41e6-a540-1f4bf92ef96c
# ╠═2ffc3d13-ab1b-4b03-864a-66a0896dbdff
# ╠═c211987a-bc3f-4ade-a5a0-e19563cf9ddd
# ╠═75641449-1801-4f22-a90b-5947f6c34b8a
# ╠═967ffab4-8abf-4447-8f09-cd0e5d624767
# ╠═01ff1c3e-6eaa-4863-bb1d-7f6f7037fafe
# ╠═d3762b94-02d7-466c-be63-a534d6a75894
# ╠═1414f977-cb37-4592-a104-217822a42ae2
# ╠═c0213d48-97da-47e3-a32c-9bf87c465c09
# ╠═07f12515-a320-434c-b1e3-73a305491f94
# ╠═e5e63932-1f89-4858-ae76-72fefa2d0139
# ╠═317602a3-d0b3-435f-af6d-112be42617d4
# ╠═a7ee49b3-3096-4818-8137-39ddc1d8195e
# ╠═86e018ba-a78c-4644-af9e-63c6fec86e40
# ╠═a51fc8b8-e2d3-4962-abd9-54a5f07d517f
# ╠═c05fca81-7f29-4576-99d6-fc35210e25b7
# ╠═d082290c-4412-42b5-aa15-e35e01923e5e
# ╠═50958e94-9de4-4b76-8d79-3b7f65b94729
