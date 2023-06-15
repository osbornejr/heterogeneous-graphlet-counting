function webpage_construction()
        @info "Building website directory structure..."
        #set up output directories
        output_dir = "$(params["website"]["website_dir"])/_assets/$(params["website"]["page_name"])"
        #run(`mkdir -p "$(params["website"]["website_dir"])/$(params["website"]["page_name"])"`)
        run(`mkdir -p "$output_dir/tableinput"`)
        run(`mkdir -p "$output_dir/plots"`)
        
        ##build run parameter table
        run_parameter_df = DataFrame(
                                     Test = params["test_name"], 
                                     Expression_cut_off = params["data_preprocessing"]["expression_cutoff"],
                                     Normalisation = params["data_preprocessing"]["norm_method"],
                                     Variance_cut_off = params["data_preprocessing"]["variance_percent"],
                                     Coexpression_measure = params["network_construction"]["coexpression"],
                                     Edge_threshold = params["network_construction"]["threshold"],
                                     Threshold_method = params["network_construction"]["threshold_method"])
        CSV.write("$output_dir/tableinput/run_parameters.csv",run_parameter_df)

        #load network (assumes is in cache)
        components,adj_matrix,network_counts,vertexlist,edgelist = get_network_construction()
        g = SimpleGraph(adj_matrix)
        #Network visualisation
        @info "Visualising network..."
        ##get largest component
        components = connected_components(g)
        largest = components[length.(components).==max(length.(components)...)]
        adj_matrix_comp = adj_matrix[largest[1],largest[1]]
        g_comp = Graph(adj_matrix_comp)
        ##update vertexlist
        vertexlist_comp = vertexlist[largest[1]]
        ##plot (either connected component or whole network
        nodefillc = [colorant"lightseagreen", colorant"orange"][(vertexlist_comp.=="coding").+1]
        draw(SVG("$output_dir/component_network.svg",16cm,16cm),gplot(g_comp,nodefillc = nodefillc))
        nodefillc = [colorant"lightseagreen", colorant"orange"][(vertexlist.=="coding").+1]
        draw(SVG("$output_dir/network.svg",16cm,16cm),gplot(g,nodefillc = nodefillc))
        
        #load in data
        raw_counts = get_input_data()
        clean_counts = cache_load("$(params["cache"]["clean_dir"])/clean_counts.jld2","clean counts")
        norm_counts = cache_load("$(params["cache"]["norm_dir"])/norm_counts.jld2","norm counts")
        sample_counts = cache_load("$(params["cache"]["sampling_dir"])/sample_counts.jld2","sample counts")
    

        ##plot before cut
        raw_data = data_from_dataframe(raw_counts,"data")
        DataPreprocessing.histogram(
                                    DataFrame([log2.(vec(sum(raw_data,dims=2))),raw_counts[!,:transcript_type]],
                                              [:sum,:transcript_type]),
                                    :sum,
                                    :transcript_type,
                                    "$output_dir/raw_data_histogram.svg",
                                    xaxis =" sum of expression (log2 adjusted)")

        ##plot after cut
        clean_data = data_from_dataframe(clean_counts,"data")
        DataPreprocessing.histogram(
                                    DataFrame([log2.(vec(sum(clean_data,dims=2))),clean_counts[!,:transcript_type]],
                                              [:sum,:transcript_type]),
                                    :sum,
                                    :transcript_type,
                                    "$output_dir/clean_data_cut_histogram.svg",
                                    xaxis =" sum of expression (log2 adjusted)")

        #boxplot(raw_counts,"raw_data_cleaned_boxplot.svg")
    
        ##Type breakdown table 
        ##set up csv string
        csv = "Step,Coding counts,Non-coding counts,Non-coding proportion\n"
        csv = csv*"Raw counts,"*string(size(raw_counts,1)-size(filter(:transcript_type=>x->x=="noncoding",raw_counts),1))*","*string(size(filter(:transcript_type=>x->x=="noncoding",raw_counts),1))*","*string(round(size(filter(:transcript_type=>x->x=="noncoding",raw_counts),1)/size(raw_counts,1),sigdigits=3))*"\n"
        csv = csv*"Clean counts,"*string(size(clean_counts,1)-size(filter(:transcript_type=>x->x=="noncoding",clean_counts),1))*","*string(size(filter(:transcript_type=>x->x=="noncoding",clean_counts),1))*","*string(round(size(filter(:transcript_type=>x->x=="noncoding",clean_counts),1)/size(clean_counts,1),sigdigits=3))*"\n"
        csv = csv*"Sample counts,"*string(size(sample_counts,1)-size(filter(:transcript_type=>x->x=="noncoding",sample_counts),1))*","*string(size(filter(:transcript_type=>x->x=="noncoding",sample_counts),1))*","*string(round(size(filter(:transcript_type=>x->x=="noncoding",sample_counts),1)/size(sample_counts,1),sigdigits=3))*"\n"
        csv = csv*"Network counts,"*string(size(network_counts,1)-size(filter(:transcript_type=>x->x=="noncoding",network_counts),1))*","*string(size(filter(:transcript_type=>x->x=="noncoding",network_counts),1))*","*string(round(size(filter(:transcript_type=>x->x=="noncoding",network_counts),1)/size(network_counts,1),sigdigits=3))*"\n"
        write("$output_dir/tableinput/type_representation.csv",csv)
        
        #Degrees
        #homogonous degree distribution
        degrees = vec(sum(adj_matrix,dims=2))
        


        p = plot(DataFrame([sort(degrees)],:auto),x = "x1",Geom.histogram,Guide.title("Degree distribution"),Guide.xlabel("degree"));
        draw(SVG("$output_dir/degree_distribution.svg"),p)
        
        #degrees for each transcript type
        for type in unique(vertexlist)
                p = plot(DataFrame([sort(degrees[vertexlist.==type])],:auto),x = "x1",Geom.histogram,Guide.title("Degree distribution"),Guide.xlabel("degree"));
                draw(SVG("$output_dir/$(type)_degree_distribution.svg"),p)
        end
        
        ## Hubs
        @info "Identifying hubs..."
        deg_thresh = Int(floor(mean(degrees)+2*std(degrees)))
        nodefillc = [colorant"black", colorant"red"][(degrees.>deg_thresh).+1]
        draw(SVG("$output_dir/two_std_hub_network.svg",16cm,16cm),gplot(g,nodefillc = nodefillc))
        deg_thresh = 70#mean(degrees)+2*std(degrees)
        nodefillc = [colorant"black", colorant"red"][(degrees.>deg_thresh).+1]
        draw(SVG("$output_dir/alt_hub_network.svg",16cm,16cm),gplot(g,nodefillc = nodefillc))
        
        ##Network stats table
        ##set up csv string
        csv = "Nodes,Edges,Components,Nodes in largest component,Maximal degree\n"
        csv = csv*string(length(vertexlist))*","*string(length(edgelist))*","*string(length(components))*","*string(length(largest...))*","*string(max(degrees...))*"\n"
        write("$output_dir/tableinput/network_stats.csv",csv)
end

