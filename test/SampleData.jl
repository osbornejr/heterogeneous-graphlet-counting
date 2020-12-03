
#add variance column to normalised data
variance = vec(var(data, dims=2))
insertcols!(norm_counts,"variance"=>variance)

#p = plot(x = log2.(variance), Geom.histogram(bincount = 100,density = false),Guide.xlabel("variance (log2)"),Guide.ylabel("frequency"),Guide.title("Variance histogram"));
#draw(SVG("variance.svg"),p)

#To run locally (ie on laptop) we must still cut down number of transcripts to managable level here. TODO put into function? At the moment, we select the N/m highest variable transcripts of each type, where m is the number of types. 
#N=100
#norm_counts_sample_noncoding=sort(norm_counts[norm_counts[:transcript_type].=="noncoding",:],:variance)[end-Int(N/2)+1:end,:]
#norm_counts_sample_coding=sort(norm_counts[norm_counts[:transcript_type].=="coding",:],:variance)[end-Int(N/2)+1:end,:]
#norm_counts_sample = outerjoin(norm_counts_sample_noncoding,norm_counts_sample_coding,on = names(norm_counts))
#data=Array(norm_counts_sample[!,2:13])
#

## We really do want a better method to do this, preferably one that scales easily. It has been decided to go for a filter that takes the top X% of each transcript type, to maintain the ratios of each type without being too prescriptive on what does and doesn't end up in the network.
X = 0.01
norm_counts_sample_noncoding=sort(norm_counts[norm_counts[:transcript_type].=="noncoding",:],:variance)[Int(round(end*(1-X))):end,:]
norm_counts_sample_coding=sort(norm_counts[norm_counts[:transcript_type].=="coding",:],:variance)[Int(round(end*(1-X))):end,:]
norm_counts_sample = outerjoin(norm_counts_sample_noncoding,norm_counts_sample_coding,on = names(norm_counts))
data=Array(norm_counts_sample[!,2:13])
