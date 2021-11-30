using Distributed, JLD, LightGraphs, GraphPlot, Colors, Random, Glob, Distributions

function distributed_setup(inclusions::Array{String,1})

        ##set up for distributed mode

        #first clean to make sure there are no stray workers already around
        if(length(workers())!=Threads.nthreads())
                rmprocs(workers())
                #add workers equal to the number of available cpus      
                addprocs(Threads.nthreads())
                #addprocs(8)
                @everywhere inclusions
                for inc in inclusions
                        @everywhere include(inc)
                end
        end
end
function webpage_construction(raw_counts::DataFrame,params::RunParameters)


        @info "Building directory structure..."
        ##establish output directories  
        run(`mkdir -p "$(params.website_dir)/_assets/$(params.page_name)/tableinput"`)
        run(`mkdir -p "$(params.website_dir)/_assets/$(params.page_name)/plots"`)

        run_parameter_df = DataFrame(Test = params.test_name, Expression_cut_off = params.expression_cutoff,Normalisation = params.norm_method, Variance_cut_off = params.variance_percent, Coexpression_measure = params.coexpression, Edge_threshold = params.threshold, Threshold_method = params.threshold_method)
        CSV.write("$(params.website_dir)/_assets/$(params.page_name)/tableinput/run_parameters.csv",run_parameter_df)
        #set up output directories
        run(`mkdir -p "$(params.website_dir)/$(params.page_name)"`)
        run(`mkdir -p "$(params.website_dir)/_assets/$(params.page_name)/tableinput"`)
        run(`mkdir -p "$(params.website_dir)/_assets/$(params.page_name)/plots"`)
                
        ##Cache setup
        #New method: cache directory updates folder by folder as we go. Avoids need for moving files around,symbolic links etc. 
        cache_dir = "$cwd/output/cache/$(params.test_name)"
        run(`mkdir -p $(cache_dir)`)
        ## Check if there are cached files for runs with similar parameters (i.e. the same before irelevant parameters are invoked for a specific output)
        #TODO improve this process to automatically detect via a dependency tree (will also tidy up cache dir)
        #similarity matrix
        #sim_check = glob("output/cache/$(params.test_name)_$(params.expression_cutoff)_$(params.norm_method)_$(params.variance_percent)_$(params.coexpression)*/similarity*",cwd)
        #check if any already exist
        #if(length(sim_check)>0)
                #check if actual cache already exists
        #       if (!(cache_dir in first.(splitdir.(sim_check))))
        #               run(`cp $(sim_check[1]) $(cache_dir)/`)
        #       end
        #end
        #functional annotation
        #func_check = glob("output/cache/$(params.test_name)_$(params.expression_cutoff)_$(params.norm_method)_$(params.variance_percent)_$(params.coexpression)_$(params.threshold)_$(params.threshold_method)*/simil_$(params.threshold)_$(params.threshold_method)*/func*",cwd)
        #check if any already exist
        #if(length(func_check)>0)
                #check if actual cache already exists
        #       if (!(cache_dir in first.(splitdir.(func_check))))
        #               run(`cp $(func_check[1]) $(cache_dir)/`)
        #       end
        #end
        #Processing data:
        @info "Processing raw data..."
        raw_data = Array(select(raw_counts,filter(x->occursin("data",x),names(raw_counts))))
        
        ## Clean - remove transcripts with total counts across all samples less than Cut
        ##update cache
        cache_dir = "$cwd/output/cache/$(params.test_name)/cutoff/$(params.expression_cutoff)"
        run(`mkdir -p $(cache_dir)`)
        ##plot before cut
        histogram(DataFrame([log2.(vec(sum(raw_data,dims=2))),raw_counts[!,:transcript_type]],[:sum,:transcript_type]),:sum,:transcript_type,"$(params.website_dir)/_assets/$(params.page_name)/raw_data_histogram.svg",xaxis =" sum of expression (log2 adjusted)")
        
        clean_counts=raw_counts[vec(sum(raw_data,dims = 2 ).>=params.expression_cutoff),:]
        clean_data = Array(select(clean_counts,filter(x->occursin("data",x),names(clean_counts))))
        
        ##plot after cut
        histogram(DataFrame([log2.(vec(sum(clean_data,dims=2))),clean_counts[!,:transcript_type]],[:sum,:transcript_type]),:sum,:transcript_type,"$(params.website_dir)/_assets/$(params.page_name)/clean_data_cut_histogram.svg",xaxis =" sum of expression (log2 adjusted)")
        
        #boxplot(raw_counts,"raw_data_cleaned_boxplot.svg")
        
        ### Normalisation
        ##update cache
        cache_dir = "$cwd/output/cache/$(params.test_name)/cutoff/$(params.expression_cutoff)/normalisation/$(params.norm_method)"
        run(`mkdir -p $(cache_dir)`)
        norm_data=library_size_normalisation(clean_data,params.norm_method)
        norm_counts = copy(clean_counts)
        norm_counts[:,findall(x->occursin("data",x),names(norm_counts))] = norm_data
        
        ##Sampling for most variable transcripts
        ##update cache
        cache_dir = "$cwd/output/cache/$(params.test_name)/cutoff/$(params.expression_cutoff)/normalisation/$(params.norm_method)/sampling/$(params.variance_percent)"
        run(`mkdir -p $(cache_dir)`)
        #add variance column to normalised data
        variance = vec(var(norm_data, dims=2))
        norm_counts.variance = variance
        sample_counts_noncoding=sort(norm_counts[norm_counts[:transcript_type].=="noncoding",:],:variance)[Int(round(end*(1-params.variance_percent))):end,:]
        sample_counts_coding=sort(norm_counts[norm_counts[:transcript_type].=="coding",:],:variance)[Int(round(end*(1-params.variance_percent))):end,:]
        sample_counts = outerjoin(sample_counts_noncoding,sample_counts_coding,on = names(norm_counts))
        sample_data = Array(select(sample_counts,filter(x->occursin("data",x),names(sample_counts))))
        ##Network construction
        @info "Constructing the network..."
        ##Measure of coexpression
        ##update cache
        cache_dir = "$cwd/output/cache/$(params.test_name)/cutoff/$(params.expression_cutoff)/normalisation/$(params.norm_method)/sampling/$(params.variance_percent)/similarity/$(params.coexpression)"
        run(`mkdir -p $(cache_dir)`)
        #similarity_matrix=mutual_information(data)
        ## file to cache similarity matrix for use later:
        sim_file = "$cache_dir/similarity_matrix.jld"
        if (isfile(sim_file))
                @info "Loading similarity matrix from $cache_dir..."
                similarity_matrix = JLD.load(sim_file,"similarity_matrix")
        else
                @info "Generating similarity matrix..."
                similarity_matrix = coexpression_measure(sample_data,params.coexpression)
                @info "Saving similarity matrix at $cache_dir..."
                JLD.save(sim_file,"similarity_matrix",similarity_matrix)
        end
        
        ## Adjacency matrix 
        ##update cache
        cache_dir = "$cwd/output/cache/$(params.test_name)/cutoff/$(params.expression_cutoff)/normalisation/$(params.norm_method)/sampling/$(params.variance_percent)/similarity/$(params.coexpression)/threshold/$(params.threshold)/threshold_method/$(params.threshold_method)"
        run(`mkdir -p $(cache_dir)`)
        adj_file = "$cache_dir/adjacency_matrix.jld"
        if (isfile(adj_file))
                @info "Loading adjacency matrix from $cache_dir..."
                pre_adj_matrix = JLD.load(adj_file,"pre-adj_matrix")
                adj_matrix = JLD.load(adj_file,"adjacency_matrix")
        else
                @info "Generating adjacency matrix..."
                if (params.threshold_method=="empirical_dist")
                        pre_adj_matrix = empirical_dist_adjacency(similarity_matrix,params.threshold)
                elseif (params.threshold_method=="empirical_dist_zero")
                        pre_adj_matrix = empirical_dist_zero_adjacency(similarity_matrix,params.threshold)
                elseif (params.threshold_method=="hard")
                        pre_adj_matrix = adjacency(similarity_matrix,params.threshold)
                elseif (params.threshold_method=="top")
                        ##TODO setting top x value here for now; should be a parameter, but as an Int rather than Float as threshold param is for other methods
                        pre_adj_matrix = top_adjacency(similarity_matrix,10)
                end
                ##form final adjacency matrix
                adj_matrix = copy(pre_adj_matrix)
                adj_matrix = adj_matrix[:,vec(sum(pre_adj_matrix,dims=1).!=0)]
                adj_matrix = adj_matrix[vec(sum(pre_adj_matrix,dims=2).!=0),:]
                @info "Saving adjacency matrix at $cache_dir..."
                JLD.save(adj_file,"adjacency_matrix",adj_matrix,"pre-adj_matrix",pre_adj_matrix)
        end
        #Trim nodes with degree zero
        network_counts = sample_counts[vec(sum(pre_adj_matrix,dims=2).!=0),:]
        network_data = sample_data[vec(sum(pre_adj_matrix,dims=2).!=0),:]
        #maintain list of vertices in graph
        vertex_names = network_counts[:transcript_id]
        vertexlist = copy(network_counts[:transcript_type])     
        edgelist = edgelist_from_adj(adj_matrix)
        
        #Synthetic test (just override vertex and edge lists here-- is that ok?)
        if(params.test_name == "Synthetic")
            n = length(vertexlist) 
            m = length(edgelist) 
            # construct erdos renyi random network based on vertex and edge structure of real network
            edgelist = Pair.(collect(edges(erdos_renyi(n,m/(n*(n-1)/2)))))
            #percentage of coding vertices in synthetic network
            percentage = 0.72
            vertexlist = vcat(repeat(["coding"],Int(floor(percentage*n))),repeat(["noncoding"],Int(ceil((1-percentage)*n))))
        end

        #Network visualisation
        g = SimpleGraph(adj_matrix)
        ##get largest component
        components = connected_components(g)
        largest = components[length.(components).==max(length.(components)...)]
        adj_matrix_comp = adj_matrix[largest[1],largest[1]]
        g_comp = Graph(adj_matrix_comp)
        ##update vertexlist
        vertexlist_comp = vertexlist[largest[1]]

        if (params.visualise==true)
                @info "Visualising network..."
                ##plot (either connected component or whole network
                nodefillc = [colorant"lightseagreen", colorant"orange"][(vertexlist_comp.=="coding").+1]
                draw(SVG("$(params.website_dir)/_assets/$(params.page_name)/component_network.svg",16cm,16cm),gplot(g_comp,nodefillc = nodefillc))
                nodefillc = [colorant"lightseagreen", colorant"orange"][(vertexlist.=="coding").+1]
                draw(SVG("$(params.website_dir)/_assets/$(params.page_name)/network.svg",16cm,16cm),gplot(g,nodefillc = nodefillc))
        end
        
        #Network Analysis
        ##update cache-- here we split cache into separate analysis folders, as each analysis is independent of the others
        cache_dir = "$cwd/output/cache/$(params.test_name)/cutoff/$(params.expression_cutoff)/normalisation/$(params.norm_method)/sampling/$(params.variance_percent)/similarity/$(params.coexpression)/threshold/$(params.threshold)/threshold_method/$(params.threshold_method)/analysis"
        run(`mkdir -p $(cache_dir)`)
        @info "Analysing network..."
        #Type representations 
        ##set up csv string
        csv = "Step,Coding counts,Non-coding counts,Non-coding proportion\n"
        csv = csv*"Raw counts,"*string(size(raw_counts,1)-size(filter(:transcript_type=>x->x=="noncoding",raw_counts),1))*","*string(size(filter(:transcript_type=>x->x=="noncoding",raw_counts),1))*","*string(round(size(filter(:transcript_type=>x->x=="noncoding",raw_counts),1)/size(raw_counts,1),sigdigits=3))*"\n"
        csv = csv*"Clean counts,"*string(size(clean_counts,1)-size(filter(:transcript_type=>x->x=="noncoding",clean_counts),1))*","*string(size(filter(:transcript_type=>x->x=="noncoding",clean_counts),1))*","*string(round(size(filter(:transcript_type=>x->x=="noncoding",clean_counts),1)/size(clean_counts,1),sigdigits=3))*"\n"
        csv = csv*"Sample counts,"*string(size(sample_counts,1)-size(filter(:transcript_type=>x->x=="noncoding",sample_counts),1))*","*string(size(filter(:transcript_type=>x->x=="noncoding",sample_counts),1))*","*string(round(size(filter(:transcript_type=>x->x=="noncoding",sample_counts),1)/size(sample_counts,1),sigdigits=3))*"\n"
        csv = csv*"Network counts,"*string(size(network_counts,1)-size(filter(:transcript_type=>x->x=="noncoding",network_counts),1))*","*string(size(filter(:transcript_type=>x->x=="noncoding",network_counts),1))*","*string(round(size(filter(:transcript_type=>x->x=="noncoding",network_counts),1)/size(network_counts,1),sigdigits=3))*"\n"
        write("$(params.website_dir)/_assets/$(params.page_name)/tableinput/type_representation.csv",csv)
        
        #Degrees
        #homogonous degree distribution
        degrees = vec(sum(adj_matrix,dims=2))
        p = plot(DataFrame([sort(degrees)]),x = "x1",Geom.histogram,Guide.title("Degree distribution"),Guide.xlabel("degree"));
        draw(SVG("$(params.website_dir)/_assets/$(params.page_name)/degree_distribution.svg"),p)
        
        #degrees for each transcript type
        for type in unique(vertexlist)
                p = plot(DataFrame([sort(degrees[vertexlist.==type])]),x = "x1",Geom.histogram,Guide.title("Degree distribution"),Guide.xlabel("degree"));
                draw(SVG("$(params.website_dir)/_assets/$(params.page_name)/$(type)_degree_distribution.svg"),p)
        end
        
        ## Hubs
        if (params.visualise==true)
                @info "Identifying hubs..."
                deg_thresh = Int(floor(mean(degrees)+2*std(degrees)))
                nodefillc = [colorant"black", colorant"red"][(degrees.>deg_thresh).+1]
                draw(SVG("$(params.website_dir)/_assets/$(params.page_name)/two_std_hub_network.svg",16cm,16cm),gplot(g,nodefillc = nodefillc))
                deg_thresh = 70#mean(degrees)+2*std(degrees)
                nodefillc = [colorant"black", colorant"red"][(degrees.>deg_thresh).+1]
                draw(SVG("$(params.website_dir)/_assets/$(params.page_name)/alt_hub_network.svg",16cm,16cm),gplot(g,nodefillc = nodefillc))
        end     
        ## Community structure
        ##update cache
        anal_dir = "$cache_dir/communities"
        run(`mkdir -p $(anal_dir)`)
        @info "Identifying communities..."
        ##use gene ids here, as they have more chance of getting a GO annotation
        if(params.func_annotate==true)
                vertex_gene_names = network_counts[:gene_id]
                community_vertices = get_community_structure(adj_matrix,vertex_gene_names,"louvain",threejs_plot = true,plot_prefix = "$(params.website_dir)/$(params.page_name)") 
                ## functional annotations of communities
                func_file = "$anal_dir/func_annotations.jld" 
                if (isfile(func_file))
                        functional_annotations = JLD.load(func_file,"functional annotations")
                else
                        functional_annotations = get_functional_annotations(community_vertices,ensembl_version = "75",write_csv = true, csv_dir ="$(params.website_dir)/_assets/$(params.page_name)/tableinput/")       
                        JLD.save(func_file,"functional annotations",functional_annotations)
                end
        else
                community_vertices = get_community_structure(adj_matrix,vertex_names,"louvain",threejs_plot = true,plot_prefix = "$(params.website_dir)/$(params.page_name)") 
        end
        #find nodes who are not in any community (usually because they are not in connected component).
        community_orphans = findall(x->x==0,in.(1:length(vertexlist),Ref(community_vertices.name)))
        ##Network stats table
        ##set up csv string
        csv = "Nodes,Edges,Components,Nodes in largest component,Maximal degree,communities detected\n"
        csv = csv*string(length(vertexlist))*","*string(length(edgelist))*","*string(length(components))*","*string(length(largest...))*","*string(max(degrees...))*","*string(length(unique(community_vertices.group)))*"\n"
        write("$(params.website_dir)/_assets/$(params.page_name)/tableinput/network_stats.csv",csv)

        if(params.graphlet_counting==true)

                ##update cache
                anal_dir = "$cache_dir/graphlets"
                run(`mkdir -p $(anal_dir)`)
                graphlet_file = "$anal_dir/graphlets.jld" 
                if (isfile(graphlet_file))
                        @info "Loading graphlet counts from $anal_dir..."
                        graphlet_counts = JLD.load(graphlet_file,"graphlets")
                        timer = JLD.load(graphlet_file,"time")
                else
                        @info "Counting graphlets..."
                        timer=@elapsed graphlet_counts = count_graphlets(vertexlist,edgelist,4,run_method="distributed-old")
                        #graphlet_concentrations = concentrate(graphlet_counts) 
                        @info "Saving graphlet counts at $anal_dir..."
                        ##save the per-edge array as well in case we need it in the future (exp for debugging)
                        JLD.save(graphlet_file,"graphlets",graphlet_counts,"time",timer)
        
                end
        

                #method to deduce which run method the null model should use given the run time of graphlets above ($timer)
                if (timer<30)
                        null_run = "distributed-short"
                else
                        null_run = "distributed-long"
                end

                ##Typed representations
                @info "Looking at typed representations of graphlets..."
                
                ##update cache
                rep_dir = "$anal_dir/typed_representations/nullmodel/$(params.null_model_size)_simulations"
                run(`mkdir -p $(rep_dir)`)
                N=params.null_model_size
                rand_graphlets_file = "$rep_dir/rand_graphlets.jld"
                if (isfile(rand_graphlets_file))
                        @info "Loading randomised vertices and graphlet counts from $rep_dir..."
                        rand_types_set = JLD.load(rand_graphlets_file,"rand vertices")
                        rand_graphlet_collection = JLD.load(rand_graphlets_file,"rand graphlets")
                else
                        ## randomise node types
                        
                        #number of randomised graphs
                        rand_types_set = [copy(vertexlist) for i in 1:N]
                        broadcast(shuffle!,rand_types_set) 
                        @info "Counting graphlets on null model" 
                        if (null_run=="distributed-short")
                                #rand_graphlet_counts = count_graphlets.(rand_types_set,Ref(edgelist),4,run_method="distributed-old")
                                rand_graphlet_counts = @showprogress pmap(x->count_graphlets(x,edgelist,4,run_method="serial"),rand_types_set,batch_size =10)
                        end
                        if (null_run=="distributed-long")
                                #rand_graphlet_counts = count_graphlets.(rand_types_set,Ref(edgelist),4,run_method="distributed")
                                rand_graphlet_counts = @showprogress map(x->count_graphlets(x,edgelist,4,run_method="distributed"),rand_types_set)
                        end
                        rand_graphlet_collection = vcat(collect.(rand_graphlet_counts)...)
                        @info "Saving random graphlet count information at $rep_dir..."
                        JLD.save(rand_graphlets_file,"rand graphlets",rand_graphlet_collection,"rand vertices",rand_types_set)
                end
                
                
                rand_df = DataFrame(graphlet = broadcast(first,rand_graphlet_collection),value = broadcast(last,rand_graphlet_collection))
                real_df = DataFrame(graphlet = broadcast(first,collect(graphlet_counts)),value = broadcast(last,collect(graphlet_counts)))
                
                ##function to get the n permutations of a set xs
                all_perm(xs, n) = vec(map(collect, Iterators.product(ntuple(_ -> xs, n)...)))
                
                ##convert graphlet_counts dict output to default dictionary, returning 0 for graphlets that don't exist in the real network
                real_dict = DefaultDict(0,graphlet_counts)
                hom_graphlets = unique(last.(split.(unique(real_df[:graphlet]),"_")))
                ##array to store all homogonous graphlet dfs
                hog_array = Array{DataFrame,1}(undef,length(hom_graphlets))       
                hog_array_under = Array{DataFrame,1}(undef,length(hom_graphlets)) 
                ##array to store all homgoenous summaries
                merged_summaries= Array{DataFrame,1}(undef,length(hom_graphlets))       
                for (i,hog) in enumerate(hom_graphlets)
                        hog_df= DataFrame()
                        hog_df_under= DataFrame()
                        ## restrict info to just hg 
                        real_fil = filter(:graphlet=>x->occursin(hog,x),real_df)
                        rand_fil = filter(:graphlet=>x->occursin(hog,x),rand_df)
                
                        #get hetero subgraphlets within homogonous type (problem: might not be complete set present in real/rand outputs?)
                        if (occursin("4",hog))
                                het_graphlets = union(first.(split.(real_fil[:graphlet],"_4")),first.(split.(rand_fil[:graphlet],"_4")))
                        elseif (occursin("3",hog))
                                het_graphlets = union(first.(split.(real_fil[:graphlet],"_3")),first.(split.(real_fil[:graphlet],"_3")))
                        end
                        #store summaries (for TikZ plot)
                        summaries = DataFrame()
                        for heg in het_graphlets
                                rand_fil_fil = filter(:graphlet=>x->x==heg*"_"*hog,rand_df)
                                transform!(rand_fil_fil,:value =>ByRow(x-> log(x))=>:log_value)
                                
                                ##get summary (for tikZ plot)
                                summary = describe(rand_fil_fil[:,3:3],:min,:q25,:median,:q75,:max)
                                summary.variable = heg*"_"*hog
                                append!(summaries,summary)
                                ##histogram for each heterogeneous graphlet
                                histogram(rand_fil_fil,:value,:graphlet,"$(params.website_dir)/_assets/$(params.page_name)/plots/$(heg)_$(hog)_histogram.svg")
                                ##log version
                                histogram(rand_fil_fil,:log_value,:graphlet,"$(params.website_dir)/_assets/$(params.page_name)/plots/$(heg)_$(hog)_log_histogram.svg")
                                rand_vals = rand_fil_fil[!,:value]
                                log_rand_vals =rand_fil_fil[!,:log_value]
                                rand_exp = sum(rand_vals)/N
                                log_rand_exp = sum(log_rand_vals)/N
                                real_obs = real_dict[heg*"_"*hog]
                                log_real_obs = log(real_dict[heg*"_"*hog]+1)
                                ## Monte Carlo method
                                r = sum(rand_vals.>=real_obs)
                                r_under = sum(rand_vals.<=real_obs)
                                ##special case: real_obs is zero.. Zero occurences do not register in rand_vals table, but clearly every rand network has equal or greater count 
                                if(real_obs==0)
                                        r = N
                                end
                                p_value = (r+1)/(N+1)
                                p_value_under = (r_under+1)/(N+1)                       
                                ##Z-score method (assumes either normal or lognormal distribution of counts in sims)
                                #using real values
                        #       z_score = (abs(real_obs) - rand_exp)/std(rand_vals)
                                ##using log values
                                #z_score = (abs(log_real_obs) - log_rand_exp)/std(log_rand_vals)
                                #p_value =pdf(Normal(0,1),z_score) 
                                append!(hog_df,DataFrame(Graphlet = heg*"_"*hog, Expected = rand_exp,Observed = real_obs, p_value = p_value))   
                                append!(hog_df_under,DataFrame(Graphlet = heg*"_"*hog, Expected = rand_exp,Observed = real_obs, p_value = p_value_under))       
                        end
                        ##take log values to plot
                        log_real_fil = copy(real_fil)
                        log_real_fil.value =log.(log_real_fil.value)
                        log_rand_fil = copy(rand_fil)
                        log_rand_fil.value =log.(log_rand_fil.value)
                        #SVG plot (for web)
                        p = plot(layer(filter(:graphlet=>x->occursin(hog,x),log_real_fil),x = :graphlet,y = :value, Geom.point,color=["count in graph"]),Guide.xticks(label=true),Theme(key_position = :none),Guide.xlabel(nothing),Guide.ylabel("log value"),Guide.yticks(orientation=:vertical),layer(filter(:graphlet=>x->occursin(hog,x),log_rand_fil),x=:graphlet,y=:value,Geom.boxplot(suppress_outliers = true),color=:graphlet));
                        draw(SVG("$(params.website_dir)/_assets/$(params.page_name)/plots/$(hog)_boxplot.svg",4inch,6inch),p)
                        #TeX plot (via PGFPlots) 
                        #add real log values to summaries, order from lowest to highest)
                        summaries.values = log_real_fil.value
                        sort!(summaries,:values)
                        tex_boxplot(summaries[!,Not(:values)],summaries.values,"output/share/$(hog)_boxplot.tex","input",ylabel="")
                        ##leave real values attached here, improves tikz layout
                        merged_summaries[i] = summaries
                        hog_array[i] = hog_df
                        hog_array_under[i] = hog_df_under
                end

                #merged boxplots
                tex_merged_boxplot(merged_summaries,"output/share/merged_boxplot.tex","input",ylabel = "log value")

                ## find significant graphlets
                sig_graphlets = vcat(filter.(:p_value=>x->x<0.05,hog_array)...)
                insig_graphlets = vcat(filter.(:p_value=>x->x<0.05,hog_array_under)...)
                
                #save in output cache
                html_table_maker(sig_graphlets,"$rep_dir/sig_type_representations.html",imgs=sig_graphlets.Graphlet)                          
                html_table_maker(insig_graphlets,"$rep_dir/insig_type_representations.html",imgs=sig_graphlets.Graphlet)                              
                #save for website version
                html_table_maker(sig_graphlets,"$(params.website_dir)/_assets/$(params.page_name)/sig_type_representations.html",imgs=sig_graphlets.Graphlet,figpath="../figs/")
                html_table_maker(insig_graphlets,"$(params.website_dir)/_assets/$(params.page_name)/insig_type_representations.html",imgs=insig_graphlets.Graphlet,figpath="../figs/")  
                ##look at edge types in randomised networks
                real_type_edgecounts = countmap(splat(tuple).(sort.(eachrow(hcat(map(x->vertexlist[x],first.(edgelist)),map(x->vertexlist[x],last.(edgelist)))))))
                rand_types_edgecounts = map(y->(countmap(splat(tuple).(sort.(eachrow(hcat(map(x->y[x],first.(edgelist)),map(x->y[x],last.(edgelist)))))))),rand_types_set)
                rand_edge_collection = vcat(collect.(rand_types_edgecounts)...)
                rand_edge_df = DataFrame(graphlet = broadcast(first,rand_edge_collection),value = broadcast(last,rand_edge_collection))
                random_edges = DataFrame()
                for t in unique(rand_edge_df[:graphlet])
                        rand_vals = filter(:graphlet=>x->x==t,rand_edge_df)[!,:value]
                        rand_exp = sum(rand_vals)/N
                        real_obs = real_type_edgecounts[t]
                        append!(random_edges,DataFrame(Graphlet = first(t)*"_"*last(t)*"_edge", Expected = rand_exp,Observed = real_obs))       
                end
                
                #pretty_table(random_edges,backend=:html,standalone = false)
                
                #@time motif_counts = find_motifs(edgelist,"hetero_rewire",100, typed = true, typelist = vec(vertexlist),plotfile="$cache_dir/motif_detection.svg",graphlet_size = 4)

#               ### Validation steps
#               # looking at identified significant graphlets and seeing if they check out biologically
#               # the most taxing step is to identify the graphlets that are coincident in some way to KEGG pathways. We cache these coincidents as a dataframe (using CSV instead of JLD)  
                val_dir = "$anal_dir/validation"
                run(`mkdir -p $(val_dir)`)
                coincidents_file ="$val_dir/coincidents.csv"
                if (isfile(coincidents_file))
                        @info "Loading coincidents dataframe from $val_dir..."
                        Coincidents = CSV.read(coincidents_file)
                        #because CSV converts the array columns to strings, we have to convert back (cost of using the easy/dirty CSV option!)
                        fix(g) = split(replace(replace(replace(replace(g,("["=>"")),("]"=>"")),("\""=>"")),(" "=>"")),",")
                        fix_int(g) = map(x->parse(Int,x),split(replace(replace(g,("["=>"")),("]"=>"")),","))
                        fix_bool(g) = BitArray(map(x->parse(Int,x),split(replace(replace(replace(g,("["=>"")),("]"=>"")),("Bool"=>"")),",")))
                        Coincidents.Vertices = fix_int.(Coincidents.Vertices)
                        Coincidents.Entrez = fix_int.(Coincidents.Entrez)
                        Coincidents.Ensembl = fix.(Coincidents.Ensembl)
                        Coincidents.Transcript_type = fix.(Coincidents.Transcript_type)
                        Coincidents.Inclusion = fix_bool.(Coincidents.Inclusion)

                else
                        @info "Conducting per graphlet pathway coincidence analysis..."
                        Coincidents = graphlet_coincidences(vertexlist,vertex_names,"transcripts",adj_matrix)
                        @info "Saving coincidents at $val_dir..."
                        CSV.write(coincidents_file,Coincidents)
        
                end
                

                ##Coincident analysis
                #get baseline entrez and kegg info about transcripts
                entrez_id_vector, candidates = get_KEGG_pathways(vertex_names,"transcripts")
                candidate_pathways = collect(keys(candidates))
                #find which types are excluded in general, and then only as a cause of having no Entrez id
                Coincidents.excluded = [Coincidents.Transcript_type[i][Coincidents.Inclusion[i].==0] for i in 1:size(Coincidents)[1]]
                Coincidents.excluded_Entrez = [Coincidents.Transcript_type[i][Coincidents.Entrez[i].==0] for i in 1:size(Coincidents)[1]]
                #find only those coincidents that involve non-coding transcripts
                Coincidents_noncoding = Coincidents[findall(x-> "noncoding" in x, Coincidents.Transcript_type),:]
                
                #Nonuniformity test: finds graphlets where the included nodes differ across different pathways (sampling for now for speed)
                graphlet = "4-star"
                nonuniforms = []
                for y in filter(p->last(p)>1,countmap(filter(:Hom_graphlet=> x-> x == graphlet, Coincidents[1:199000,:]).Vertices))
                    test = sum(filter(:Vertices => x-> x == first(y),filter(:Hom_graphlet => x-> x == graphlet,Coincidents))[8])/last(y)
                    if (sum(((test.>0) - (test.<1)).==0)>0)
                        push!(nonuniforms,first(y))
                    end
                end
                #types of exlusions: which transcript types are most likely to be missing from the pathway in a graphlet
                countmap([Coincidents.Transcript_type[i][Coincidents.Inclusion[i].==0] for i in 1:size(Coincidents)[1]])
                countmap([Coincidents_noncoding.Transcript_type[i][Coincidents_noncoding.Inclusion[i].==0] for i in 1:size(Coincidents_noncoding)[1]])
                
                #agreement between entrez and inclusion info
                sum(map(x->x.!==0,Coincidents.Entrez).==Coincidents.Inclusion)
                sum(map(x->x.!==0,Coincidents_noncoding.Entrez).==Coincidents_noncoding.Inclusion)

                #orbit statistics
                @info "Getting orbit significance statistics..."
                #peripheral: degree one orbit
                #central: degree two orbit
                #supercentral: degree three orbit
                orbit_templates = Dict("3-path" => Dict(("central" => [0,1,0]), ("peripheral" => [1,0,1])),
                    "3-tri" => Dict(("central" => [1,1,1])),
                    "4-path" => Dict(("peripheral" => [1,0,0,1]), ("central" => [0,1,1,0])),
                    "4-star" => Dict(("peripheral" => [1,1,0,1]), ("supercentral" => [0,0,1,0])),
                    "4-tail" => Dict(("peripheral" => [0,0,0,1]), ("central" => [1,1,0,0]), ("supercentral" => [0,0,1,0])),
                    "4-cycle" => Dict(("central"=>[1,1,1,1])),
                    "4-chord" => Dict(("supercentral"=>[0,1,1,0]),("central"=>[1,0,0,1])),
                    "4-clique" => Dict(("supercentral"=>[1,1,1,1]))) 
               
                #collect dataframe for each node in this array
                #per_node_significance = Array{DataFrame,1}(undef,length(vertexlist))
                #choose just one order of graphlets (3 or 4)
                sub_Coincidents = filter(:Hom_graphlet=>x->occursin("4-",x),Coincidents)
               
                ## select whether we are looking at "detailed" or "collated" significance
                orbit_sigs_method ="collated"
                #table to showing whether each node (row) is included in each pathway (column)
                inkey = hcat([ in.(1:length(vertexlist),Ref(candidates[p])) for p in keys(candidates) ]...)
                if(orbit_sigs_method == "collated")
                    orbit_sigs = @showprogress map(x->pernode_significance(x,sub_Coincidents,collect(keys(candidates)),inkey[x,:]),1:length(vertexlist))
                end
                if(orbit_sigs_method == "detailed")
                    orbit_sigs = @showprogress map(x->pernode_significance_detail(x,sub_Coincidents,collect(keys(candidates)),inkey[x,:]),1:length(vertexlist))
                end

                ## now compare the significance profile of those nodes that are not attached to a pathway to the average pathway profile of known pathway nodes
                ##convert to array form for comparisons
                orbit_sigs_array = map(x->Array(x[2:end]),orbit_sigs)
                # this table stores the average for each pathway (of nodes that are known to be in the pathway)
                significance_bars = zeros(length(keys(candidates)),size(orbit_sigs_array[1])[2])
                for (i,c) in enumerate(keys(candidates))
                    significance_bars[i,:] = (sum(map(x->Array(x[2:end]),orbit_sigs[candidates[c]]))./length(vertexlist))[i,:]
                end
                ##now compare bars against profiles of non-pathway nodes
                ## choose a subset of nodes to look at. Can be boolean BitArray (with length equal to all nodes) or a specific list of nodes 
                subset = entrez_id_vector.==0
                putative_pathways = Array{Array{String,1}}(undef,length(orbit_sigs_array[subset]))
                for (i,t) in enumerate(orbit_sigs_array[subset])
                    putative_pathways[i] = collect(keys(candidates))[vec((sum(t.>significance_bars,dims=2).>2))]
                end
                
                # Gadfly beeswarm visualisation:
                output_dir = "output/plots/orbit_significance_$(orbit_sigs_method)/"
                run(`mkdir -p $output_dir`)
                # get data into wide format
                # size of each df
                last_col = size(orbit_sigs[1])[2]-1
                wide_orbit_sigs = vcat(map(x->stack(x,2:last_col+1),orbit_sigs)...)
                #for each pathway, we map the three orbit categories side by side
                palette = ["#db63c5","#bababa","#32a852"]
                for p in keys(candidates)
                    @info "Drawing beeswarm for $p..."    
                    p_df = filter(:Pathway=>x->x == p,wide_orbit_sigs)
                    #insertcols!(p_df,:log_value =>log.(p_df.value))
                    #define colors by whether entrez id of node is in pathway
                    in_pathway = vcat(collect(eachrow(repeat(in.(1:length(vertexlist),Ref(candidates[p])),1,last_col)))...)
                    non_entrez = vcat(collect(eachrow(repeat(entrez_id_vector.==0,1,last_col)))...)
                    coloring = CategoricalArray(in_pathway-non_entrez)
                  coloring = recode(coloring,-1=>"unidentified",0=>"not in pathway",1=>"in pathway")
                    insertcols!(p_df,:color =>coloring)
                    #remove non pathway nodes
                    #filter!(:color=>x->x!="not in pathway",p_df)
                    #remove zero nodes
                    filter!(:value=>x->x!=0,p_df)
                    pl = plot(p_df,x = :variable,y = :value, color = :color,Guide.title(p),Geom.beeswarm(padding = 1mm),Theme(bar_spacing=1mm,point_size=0.5mm),Scale.color_discrete_manual(palette...));
                    draw(SVG("$(output_dir)/$(p)_beeswarm.svg",30cm,20cm),pl)
                end
                #@time motif_counts = find_motifs(edgelist,"hetero_rewire",100, typed = true, typelist = vec(vertexlist),plotfile="$cache_dir/motif_detection.svg",graphlet_size = 4)

                #High zero exploration
                ## over all transcripts 
                total_zero_proportion = sum(map(x->x.==0,orbit_sigs_array))./length(orbit_sigs_array)
                ##over a subset
                subset = candidates[candidate_pathways[1]] 
                subset_zero_proportion = sum(map(x->x.==0,orbit_sigs_array[subset]))./length(orbit_sigs_array[subset])
                #compare
                zero_comparison = subset_zero_proportion.<total_zero_proportion
                #compare for all pathways,speficically for that pathway
                zero_scores = zeros(Int,length(candidate_pathways),last_col) 
                for (i,p) in enumerate(candidate_pathways)
                    subset = candidates[p] 
                    subset_zero_proportion = sum(map(x->x.==0,orbit_sigs_array[subset]))./length(orbit_sigs_array[subset])
                    #compare
                    zero_scores[i,:] = (subset_zero_proportion.<total_zero_proportion)[i,:]
                end
                #find those pathways with majority of zero proportions below total proportions 
                zero_passes = vec(sum(zero_scores,dims=2).>(last_col/2))
                zero_candidate_pathways = candidate_pathways[zero_passes]
                zero_orbit_sigs = map(x->filter(:Pathway=>y->y in zero_candidate_pathways,x),orbit_sigs)
                zero_orbit_sigs_array = map(x->Array(x[2:end]),zero_orbit_sigs)
                zero_candidates = Dict(Pair.(zero_candidate_pathways,[candidates[x] for x in zero_candidate_pathways]))

                #uniqueness of pathway contributors (concerns of too much overlap) 
                sig_pathway_occurences = countmap(vcat([candidates[x] for x in zero_candidate_pathways]...))
                m = max(collect(values(sig_pathway_occurences))...) 
                
                m = 8
                supersharers = first.(filter(x->last(x)==m,collect(sig_pathway_occurences)))
                #for these supersharers, find the set of pathways they are involved in
                supersharer_pathways = Array{Array{String,1},1}(undef,length(supersharers))
                for (i,s) in enumerate(supersharers)
                    list = []
                    for c in zero_candidates
                        if (s in last(c))
                            push!(list,first(c))
                        end
                    end
                    supersharer_pathways[i] = list
                end
                in_group = collect(keys(countmap(vcat(supersharer_pathways...))))
                not_in_group = zero_candidate_pathways[(!).(in.(zero_candidate_pathways,Ref(collect(keys(countmap(vcat(supersharer_pathways...)))))))]
                countmap(supersharer_pathways)




                #ecdfs
                #store ecdf functions in table
                ecdf_table = Array{ECDF,2}(undef,length(keys(candidates)),last_col)
                for (i,p) in enumerate((keys(candidates)))
                    for j in 1:last_col
                        ecdf_table[i,j] = ecdf(map(x->x[i,j],orbit_sigs_array))
                    end
                end
                #first check candidate probabilities across all categories
                known_pathway_probs = Array{Array{Float64,2},1}(undef,length(keys(candidates)))
                for (i,c) in enumerate(keys(candidates))
                    known_pathway_probs[i] = hcat([map(ecdf_table[i,j],map(x->x[i,j],orbit_sigs_array[candidates[c]])) for j in 1:last_col]...)
                end
                    

                #shape plots into a grid
                ncols = 6
                dims = fldmod(length(keys(candidates)),ncols)
                plots = Array{Union{Plot,Context},2}(undef,dims[1]+(dims[2]>0),ncols)
                #plots = Array{Union{Plot,Context},1}(undef,length(keys(candidates)))
                for (i,p) in enumerate(keys(candidates))
                    #find max over all measures for pathway
                    m = max([map(x->x[i,1],orbit_sigs_array)...,map(x->x[i,2],orbit_sigs_array)..., map(x->x[i,3],orbit_sigs_array)...]...)
                    plots[i] = plot([layer(x->ecdf(map(x->x[i,j],orbit_sigs_array))(x),0,m,color=[j]) for j in 1:last_col]...,
                              Scale.color_discrete_manual("orange", "green", "purple"),
                              Guide.title(p),
                              Guide.xlabel("count"),
                              Theme(major_label_font_size=4pt,key_position=:none));
                              #Guide.colorkey(title="orbit position"),
                              #Guide.title(p));
                end
               #Add legend pane
               legend = plot(wide_orbit_sigs,color=:variable,
                             Geom.blank,
                             Scale.color_discrete_manual("orange", "green", "purple"));
                #append blank gridspots if necessary
                for i in 1:length(plots)
                    if(!isassigned(plots,i))
                        plots[i] = context()
                    end
                end
                #TODO this needs to be more generalised for any dimension of gridstack/number of pathways 
                plots[2,6] = legend;
                draw(SVG("$(output_dir)/ecdfs.svg",30cm,20cm),gridstack(plots))

        end
end
