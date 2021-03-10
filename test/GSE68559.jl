### Include all source files TODO make this occur more fluently and automatically by creating a package, and using Revise
for src in filter(x->endswith(x,".jl"),readdir("src"))
	include(ENV["JULIA_PROJECT"]*"/src/"*src)
end
samples = CSV.read.(filter(x->occursin(".txt",x),readdir("data/GSE68559_RAW",join=true)))
sample_names = replace.(filter(x->occursin(".txt",x),readdir("data/GSE68559_RAW")),"_isoforms_expr.txt"=>"").*" data"
transcript_names=select(samples[1],1)
raw_counts=rename!(hcat(transcript_names,select.(samples,Symbol("FPKM"))...,makeunique=true),["transcript_id";sample_names])

raw_counts.transcript_type = replace(x-> occursin("lnc",x) ? "noncoding" : "coding",raw_counts.transcript_id)
data = Array(select(raw_counts,filter(x->occursin("data",x),names(raw_counts))))


## Clean - remove transcripts with total counts across all samples less than X
X = 25
raw_counts=raw_counts[vec(sum(data,dims = 2 ).>=X),:]
data=Array(raw_counts[!,2:13]);

boxplot(raw_counts,"raw_data_cleaned_boxplot.svg")

### Normalisation
data=library_size_normalisation(data,"upperquartile")
#PCs,D_1,per_sample=pca(data')
#p=pca_plot(PCs,3);
#draw(SVG("output/Construction/test_per_sample.svg"),p)

## update data to normalised version
#data=per_sample
norm_counts=copy(raw_counts)
norm_counts[!,2:13]=data
