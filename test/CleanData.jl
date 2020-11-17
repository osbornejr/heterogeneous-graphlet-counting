
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

#boxplot(norm_counts,"norm_data_boxplot.svg")

#p = plot(x = log2.(vec(data)), Geom.histogram(bincount = 100,density = false),Guide.xlabel("normalised expression (log2)"),Guide.ylabel("frequency"),Guide.title("Normalised expression histogram"));
#draw(SVG("norm_expression.svg"),p)

#norm_means = sum(data,dims=2)./12

#p = plot(x = log2.(vec(norm_means)), Geom.histogram(bincount = 100,density = false),Guide.xlabel("normalised mean expression (log2)"),Guide.ylabel("frequency"),Guide.title("Normalised mean expression histogram"));
#draw(SVG("norm_mean_expression.svg"),p)


