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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CommonMark = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DataPreprocessing = "0c67aaa8-d5ff-4929-99a0-75b09377fbc9"
DataStructures = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
Gadfly = "c91e804a-d5a3-530f-b6f0-dfbca275c004"
GraphletAnalysis = "32f39a16-8143-4a50-a7e7-080c0e917f42"
GraphletCounting = "7ac45bc0-02f1-46da-ad35-65e91b15b4e1"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
NetworkConstruction = "6c2e41d2-72ae-425a-84e9-b8f08a301efb"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Portinari = "72ee7ba2-92b2-4971-a97d-28f521fe8910"
ProjectFunctions = "a8586eae-54f0-4952-9436-ba92c8ab3181"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
CommonMark = "~0.9.1"
DataFrames = "~1.7.0"
DataStructures = "~0.18.22"
Gadfly = "~1.4.1"
HypertextLiteral = "~0.9.5"
PlutoUI = "~0.7.65"
Portinari = "~0.1.2"
Revise = "~3.6.6"
StatsBase = "~0.34.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "b8eea1e39d91f4e8d7ead0d07891ab6dda648bf9"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "f7817e2e585aa6d924fd714df1e2a84be7896c60"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.3.0"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "deddd8725e5e1cc49ee205a1964256043720a6c3"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.15"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.CategoricalArrays]]
deps = ["DataAPI", "Future", "Missings", "Printf", "Requires", "Statistics", "Unicode"]
git-tree-sha1 = "1568b28f91293458345dabba6a5ea3f183250a61"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.10.8"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "1713c74e00545bfe14605d2a2be1712de8fbcb58"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.1"

[[deps.ChangesOfVariables]]
deps = ["InverseFunctions", "LinearAlgebra", "Test"]
git-tree-sha1 = "3aa4bf1532aa2e14e0374c4fd72bed9a9d0d0f6c"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.10"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "062c5e1a5bf6ada13db96a4ae4749a4c2234f521"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.9"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "403f2d8e209681fcbd9468a8514efff3ea08452e"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.29.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.CommonMark]]
deps = ["PrecompileTools"]
git-tree-sha1 = "351d6f4eaf273b753001b2de4dffb8279b100769"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.9.1"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.Compose]]
deps = ["Base64", "Colors", "DataStructures", "Dates", "IterTools", "JSON", "LinearAlgebra", "Measures", "Printf", "Random", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "b137aa32bfe5b89996f8f87825b64ac41b9f2e16"
uuid = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
version = "0.9.6"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "b19db3927f0db4151cb86d073689f2428e524576"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.10.2"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CoupledFields]]
deps = ["LazyGrids", "LinearAlgebra", "Statistics", "StatsBase"]
git-tree-sha1 = "e2f5caf2e315453d24f0c731eb9e17fa5e81589d"
uuid = "7ad07ef1-bdf2-5661-9d2b-286fd4296dac"
version = "0.3.0"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "fb61b4812c49343d7ef0b533ba982c46021938a6"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.7.0"

[[deps.DataPreprocessing]]
deps = ["CSV", "Cairo", "Compose", "DataFrames", "Gadfly", "LinearAlgebra", "RCall", "Statistics"]
path = "/home/osbornejr/app/packages/DataPreprocessing"
uuid = "0c67aaa8-d5ff-4929-99a0-75b09377fbc9"
version = "0.1.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Deno_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "82adce9b61f01f036934b6fb2b100636df941b37"
uuid = "04572ae6-984a-583e-9378-9577a1c2574d"
version = "2.2.6+0"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.Discretizers]]
deps = ["DataStructures", "SpecialFunctions", "Statistics", "StatsBase"]
git-tree-sha1 = "3ac8e0ed0f999053c897bc5d2adf1c01d16afc00"
uuid = "6e83dbb3-75ca-525b-8ae2-3751f0dd50b4"
version = "3.2.4"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["AliasTables", "ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "0b4190661e8a4e51a842070e7dd4fae440ddb7f4"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.118"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d55dffd9ae73ff72f1c0482454dcf2ec6c6c4a63"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.5+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "53ebe7511fa11d33bec688a9178fac4e49eeee00"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "797762812ed063b9b94f6cc7742bc8883bb5e69e"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.9.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "b66970a70db13f45b7e57fbda1736e1cf72174ea"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.17.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Test"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "PDMats", "SparseArrays", "Statistics"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "301b5d5d731a0654825f1f2e906990f7141a106b"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.16.0+0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.Gadfly]]
deps = ["Base64", "CategoricalArrays", "Colors", "Compose", "Contour", "CoupledFields", "DataAPI", "DataStructures", "Dates", "Distributions", "DocStringExtensions", "Hexagons", "IndirectArrays", "IterTools", "JSON", "KernelDensity", "LinearAlgebra", "Loess", "Measures", "Printf", "REPL", "Random", "Requires", "Showoff", "Statistics"]
git-tree-sha1 = "ab56d1feb312d7202015905e43ebd1ae77cb393f"
uuid = "c91e804a-d5a3-530f-b6f0-dfbca275c004"
version = "1.4.1"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "fee60557e4f19d0fe5cd169211fdda80e494f4e8"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.84.0+0"

[[deps.GraphPlot]]
deps = ["ArnoldiMethod", "Colors", "Compose", "DelimitedFiles", "Graphs", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "4a616fbb4f4df5c51447192eba327da342a329b0"
uuid = "a2cc645c-3eea-5389-862e-a155d0052231"
version = "0.6.1"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.GraphletAnalysis]]
deps = ["DataFrames", "DataPreprocessing", "GraphletCounting", "Graphs", "NetworkConstruction", "RCall"]
path = "/home/osbornejr/app/packages/GraphletAnalysis"
uuid = "32f39a16-8143-4a50-a7e7-080c0e917f42"
version = "0.1.0"

[[deps.GraphletCounting]]
deps = ["CSV", "DataFrames", "DataStructures", "Distributed", "Graphs", "LinearAlgebra", "ProgressMeter", "StatsBase"]
path = "/home/osbornejr/app/packages/GraphletCounting"
uuid = "7ac45bc0-02f1-46da-ad35-65e91b15b4e1"
version = "0.1.0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "3169fd3440a02f35e549728b0890904cfd4ae58a"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.12.1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.Hexagons]]
deps = ["Test"]
git-tree-sha1 = "de4a6f9e7c4710ced6838ca906f81905f7385fd6"
uuid = "a1b4810d-1bce-5fbd-ac56-80944d57a21f"
version = "0.2.0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InformationMeasures]]
deps = ["Discretizers"]
git-tree-sha1 = "874d48f2026e8faf3fd55c86973fd028b02cd1a0"
uuid = "96684042-fbdc-5399-9b8e-d34e539a126c"
version = "0.3.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "8594fac023c5ce1ef78260f24d1ad18b4327b420"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.4"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "0f14a5456bdc6b9731a5682f439a672750a09e48"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.0.4+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"

[[deps.InverseFunctions]]
deps = ["Dates", "Test"]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "PrecompileTools", "Requires", "TranscodingStreams"]
git-tree-sha1 = "1059c071429b4753c0c869b75c859c44ba09a526"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.5.12"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eac1206917768cb54957c65a615460d87b455fc1"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.1+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "fc8504eca188aaae4345649ca6105806bc584b70"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.37"

[[deps.Juno]]
deps = ["Base64", "Logging", "Media", "Profile"]
git-tree-sha1 = "07cb43290a840908a771552911a6274bc6c072c7"
uuid = "e5e0dc1b-0480-54bc-9374-aad01c23163d"
version = "0.8.4"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "7d703202e65efa1369de1279c162b915e245eed1"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.9"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyGrids]]
deps = ["Statistics"]
git-tree-sha1 = "f43d10fea7e448a60e92976bbd8bfbca7a6e5d09"
uuid = "7031d0ef-c40d-4431-b2f8-61a8d2f650db"
version = "0.5.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a31572773ac1b745e0343fe5e2c8ddda7a37e997"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.0+0"

[[deps.Librsvg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pango_jll", "Pkg", "gdk_pixbuf_jll"]
git-tree-sha1 = "ae0923dab7324e6bc980834f709c4cd83dd797ed"
uuid = "925c91fb-5dd6-59dd-8e8c-345e74382d89"
version = "2.54.5+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "4ab7581296671007fc33f07a721631b8855f4b1d"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.1+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "321ccef73a96ba828cd51f2ab5b9f917fa73945a"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Loess]]
deps = ["Distances", "LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "f749e7351f120b3566e5923fefdf8e52ba5ec7f9"
uuid = "4345ca2d-374a-55d4-8d30-97f9976e7612"
version = "0.6.4"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "260dc274c1bc2cb839e758588c63d9c8b5e639d1"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "3.0.5"

[[deps.Luxor]]
deps = ["Base64", "Cairo", "Colors", "DataStructures", "Dates", "FFMPEG", "FileIO", "Juno", "LaTeXStrings", "PrecompileTools", "Random", "Requires", "Rsvg"]
git-tree-sha1 = "aa3eb624552373a6204c19b00e95ce62ea932d32"
uuid = "ae8d54c2-7ccd-5906-9d76-62fc9837b5bc"
version = "3.8.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "5de60bc6cb3899cd318d80d627560fae2e2d99ae"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.0.1+1"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Media]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "75a54abd10709c01f1b86b84ec225d26e840ed58"
uuid = "e89f7d12-3494-54d1-8411-f7d8b9ae1f27"
version = "0.5.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "030ea22804ef91648f29b7ad3fc15fa49d0e6e71"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.3"

[[deps.NetworkConstruction]]
deps = ["ColorSchemes", "Colors", "DataFrames", "DataPreprocessing", "Distributed", "Graphs", "InformationMeasures", "LinearAlgebra", "Luxor", "Printf", "ProgressMeter", "RCall", "SharedArrays", "Statistics", "StatsBase"]
path = "/home/osbornejr/app/packages/NetworkConstruction"
uuid = "6c2e41d2-72ae-425a-84e9-b8f08a301efb"
version = "0.1.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9216a80ff3682833ac4b733caa8c00390620ba5d"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.0+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "275a9a6d85dc86c24d03d1837a0010226a96f540"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.56.3+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlutoDevMacros]]
deps = ["HypertextLiteral", "InteractiveUtils", "MacroTools", "Markdown", "Random", "Requires"]
git-tree-sha1 = "b4b23b981704ac3e2c771a389c2899e69306c091"
uuid = "a0499f29-c39b-4c5c-807c-88074221b949"
version = "0.4.8"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3151a0c8061cc3f887019beebf359e6c4b3daa08"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.65"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.Portinari]]
deps = ["AbstractPlutoDingetjes", "Deno_jll", "HypertextLiteral", "InteractiveUtils", "Markdown", "Parameters", "PlutoDevMacros"]
git-tree-sha1 = "f2c25f678001ea2e0063786d2b183d9410826ed9"
uuid = "72ee7ba2-92b2-4971-a97d-28f521fe8910"
version = "0.1.2"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "66b20dd35966a748321d3b2537c4584cf40387c7"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.3.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "13c5103482a8ed1536a54c08d0e742ae3dca2d42"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.4"

[[deps.ProjectFunctions]]
deps = ["CSV", "CategoricalArrays", "Colors", "Compose", "DataFrames", "DataPreprocessing", "DataStructures", "Dates", "Distributed", "Gadfly", "GraphPlot", "GraphletAnalysis", "GraphletCounting", "Graphs", "JLD2", "NetworkConstruction", "Pkg", "ProgressMeter", "Random", "StatsBase", "YAML"]
path = "/home/osbornejr/app/packages/ProjectFunctions"
uuid = "a8586eae-54f0-4952-9436-ba92c8ab3181"
version = "0.1.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

[[deps.RCall]]
deps = ["CategoricalArrays", "Conda", "DataFrames", "DataStructures", "Dates", "Libdl", "Preferences", "REPL", "Random", "Requires", "StatsModels", "WinReg"]
git-tree-sha1 = "db17ec90d9f904b79e7877a764fdf95ff5c5f315"
uuid = "6f49c342-dc21-5d91-9882-a32aef131414"
version = "0.14.6"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "c4171a96893f861e4872b9868f5a7b19da8e9bc8"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.6.6"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.Rsvg]]
deps = ["Cairo", "Glib_jll", "Librsvg_jll"]
git-tree-sha1 = "3d3dc66eb46568fb3a5259034bfc752a0eb0c686"
uuid = "c4c386cf-5103-5370-be45-f3a111cca3b8"
version = "1.0.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "712fb0231ee6f9120e005ccd56297abbc053e7e0"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.8"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.ShiftedArrays]]
git-tree-sha1 = "503688b59397b3307443af35cd953a13e8005c16"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "2.0.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "41852b8679f78c8d8961eeadc8f62cef861a52e3"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "0feb6b9031bd5c51f9072393eb5ab3efd31bf9e4"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.13"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "29321314c920c26684834965ec2ce0dacc9cf8e5"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.4"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "35b09e80be285516e52c9054792c884b9216ae3c"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.4.0"

[[deps.StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "REPL", "ShiftedArrays", "SparseArrays", "StatsAPI", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "9022bcaa2fc1d484f1326eaa4db8db543ca8c66d"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.7.4"

[[deps.StringEncodings]]
deps = ["Libiconv_jll"]
git-tree-sha1 = "b765e46ba27ecf6b44faf70df40c57aa3a547dcb"
uuid = "69024149-9ee7-55f6-a4c4-859efe599b68"
version = "0.3.7"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.URIs]]
git-tree-sha1 = "24c1c558881564e2217dcf7840a8b2e10caeb0f9"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WinReg]]
git-tree-sha1 = "cd910906b099402bcc50b3eafa9634244e5ec83b"
uuid = "1b915085-20d7-51cf-bf83-8f477d6f5128"
version = "1.0.0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "b8b243e47228b4a3877f1dd6aee0c5d56db7fcf4"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.6+1"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fee71455b0aaa3440dfdd54a9a36ccef829be7d4"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.1+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.YAML]]
deps = ["Base64", "Dates", "Printf", "StringEncodings"]
git-tree-sha1 = "2f58ac39f64b41fb812340347525be3b590cce3b"
uuid = "ddb6d928-2868-570f-bddf-ab3f9cf99eb6"
version = "0.4.14"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.gdk_pixbuf_jll]]
deps = ["Artifacts", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Xorg_libX11_jll", "libpng_jll"]
git-tree-sha1 = "895f21b699121d1a57ecac57e65a852caf569254"
uuid = "da03df04-f53b-5353-a52f-6a8b0620ced0"
version = "2.42.13+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "522c1df09d05a71785765d19c9524661234738e9"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.11.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "cd155272a3738da6db765745b89e466fa64d0830"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.49+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d5a767a3bb77135a99e433afe0eb14cd7f6914c3"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"
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
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
