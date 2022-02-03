#using Distributed, JLD2, CSV
#using StatsBase,Gadfly,Compose,DataFrames,YAML,Dates, LightGraphs, Colors, Random, Distributions, ProgressMeter
using Pkg,ProgressMeter,DataFrames,YAML,Distributed,JLD2,CSV,StatsBase,Random,LightGraphs,Dates,Colors,Gadfly,Compose,DataStructures,CategoricalArrays



function run_all(config_file::String)
   
    if (length(workers())!=Threads.nthreads()) 
        Pkg.resolve()
        Pkg.precompile()
        @info "Setting up worker processes"
        distributed_setup(:ProjectFunctions,:GraphletCounting,:GraphletAnalysis,:NetworkConstruction)
    end
    @info "Loading parameters"
    load_config(config_file)
    @info "Loading raw counts"
    raw_counts = get_input_data()
    @info "Preprocessing raw counts"
    processed_counts = data_preprocessing(raw_counts)
    @info "Constructing network"
    adj_matrix,network_counts,vertexlist,edgelist = network_construction(processed_counts)

    cache_update("analysis")
    # @info "Finding communities"
    # com_anal = community_analysis(network_counts,adj_matrix)
    @info "Counting graphlets"
    graphlet_counts,timer = graphlet_counting(vertexlist,edgelist)

    cache_update("graphlets")
    # @info "Comparing typed graphlet representations"
    #typed_anal = typed_representations(graphlet_counts,timer,vertexlist,edgelist)
    @info "Conducting coincident graphlet analysis"
    coinc_anal = coincident_analysis(adj_matrix,network_counts,vertexlist,edgelist)
end

function load_config(config_file::String)
        #setup parameters to be read.writable for all functions in this module
        global params = YAML.load_file(config_file)
        cache_setup()
end

function cache_setup()
    ##create directory for run
    cache_dir = join([cwd,params["cache"]["base_dir"],params["test_name"]],"/") 
    run(`mkdir -p $(cache_dir)`)
    params["cache"]["cur_dir"] = cache_dir
    ##check if whole cache_dir should be removed
    if(params["cache"]["clear"]["all"])
        if(params["cache"]["clear"]["archive"])
            #if both clear cache and archive are true, then cache will be moved to archive. otherwise, cache will be deleted permanently.
            @info "Archiving previous cache files at archives/$(params["test_name"])..."
            run(`mkdir -p $(join([cwd,params["cache"]["base_dir"]],"/"))/archive`)
            run(`mv $(cache_dir) $(join([cwd,params["cache"]["base_dir"]],"/"))/archive/$(params["test_name"]*"_"*string(now()))`)
            run(`mkdir -p $(cache_dir)`)
        else
            @info "Clearing previous cache..."
            run(`rm -r $(cache_dir)`)
            run(`mkdir -p $(cache_dir)`)
        end
    end

end

function cache_update(general::String,specific::Any)
    return cache_update(general,string(specific))
end

function cache_update(general::String,specific::String="",side_dir::String="")
    ##method to update or add cache path 
    
    ##update current dir, possibly with specific subdir
    if (side_dir=="") 
        cache_dir = params["cache"]["cur_dir"]*"/"*general  
        params["cache"]["cur_dir"] = cache_dir
        run(`mkdir -p $(cache_dir)`)
        if(specific!="")
            cache_dir = params["cache"]["cur_dir"]*"/"*specific
            params["cache"]["cur_dir"] = cache_dir 
            run(`mkdir -p $(cache_dir)`)
        end
    else
        ##create side dir, possiblly with specific subdir
        cache_dir = params["cache"]["cur_dir"]*"/"*general  
        params["cache"][side_dir] = cache_dir 
        run(`mkdir -p $(cache_dir)`)
        if(specific!="")
            cache_dir = params["cache"][sid_dir]*"/"*specific
            params["cache"][side_dir] = cache_dir
            run(`mkdir -p $(cache_dir)`)
        end
    end
end

function get_input_data()
    file = join([params["cache"]["cur_dir"],"raw_counts.jld2"],"/")
    #cache_check(file)
    if isfile(file) 
        raw_counts = cache_load(file,"raw counts")  
    else 
        ##run external julia script provided in paramaters to create raw_counts cache. This is a bit messy as the script must match the file location here (and JLD2 file name)
        #but is best solution for now.
        run(`julia --project=$cwd/Project.toml $(params["raw_counts_script"])`)
        raw_counts = cache_load(file,"raw counts")  
     end
    return raw_counts  
end 

function data_preprocessing(raw_counts::DataFrame)
    
    ## Clean - remove transcripts with total counts across all samples less than Cut
    cache_update("expression_cutoff",string(params["data_preprocessing"]["expression_cutoff"]))
    file = "$(params["cache"]["cur_dir"])/clean_counts.jld2"

    ##plot before cut
    #DataPreprocessing.histogram(DataFrame([log2.(vec(sum(raw_data,dims=2))),raw_counts[!,:transcript_type]],[:sum,:transcript_type]),:sum,:transcript_type,"$(params["website"]["website_dir"])/_assets/$(params["website"]["page_name"])/raw_data_histogram.svg",xaxis =" sum of expression (log2 adjusted)")
    if(isfile(file))
        clean_counts = cache_load(file,"clean counts")
    else

        clean_counts = DataPreprocessing.clean_raw_counts(raw_counts,params["data_preprocessing"]["expression_cutoff"])
        @info "Saving clean counts to $file"
        cache_save(file,"clean counts"=>clean_counts)
    end

    ##plot after cut
    #DataPreprocessing.histogram(DataFrame([log2.(vec(sum(clean_data,dims=2))),clean_counts[!,:transcript_type]],[:sum,:transcript_type]),:sum,:transcript_type,"$(params["website"]["website_dir"])/_assets/$(params["website"]["page_name"])/clean_data_cut_histogram.svg",xaxis =" sum of expression (log2 adjusted)")

    #boxplot(raw_counts,"raw_data_cleaned_boxplot.svg")

    ### Normalisation
    ##update cache
    cache_update("normalisation",params["data_preprocessing"]["norm_method"])
    file = "$(params["cache"]["cur_dir"])/norm_counts.jld2"

    if(isfile(file))
        norm_counts = cache_load(file,"norm counts")
    else
        @info "Saving norm counts to $file"
        norm_counts = DataPreprocessing.normalise_clean_counts(clean_counts,params["data_preprocessing"]["norm_method"])
        cache_save(file,"norm counts"=>norm_counts)
    end

    ##Sampling for most variable transcripts
    ##update cache
    cache_update("sampling",string(params["data_preprocessing"]["variance_percent"]))
    file = "$(params["cache"]["cur_dir"])/sample_counts.jld2"

    if(isfile(file))
        sample_counts = cache_load(file,"sample counts")
    else
        @info "Saving sample counts to $file"
        sample_counts = DataPreprocessing.sample_norm_counts(norm_counts,params["data_preprocessing"]["variance_percent"])
        cache_save(file,"sample counts"=>sample_counts)
    end

    return sample_counts
end 

function  network_construction(sample_counts::DataFrame)

    ##Network construction
    ##Measure of coexpression
    ##update cache

    cache_update("similarity",params["network_construction"]["coexpression"])

    #similarity_matrix=mutual_information(data)
    ## file to cache similarity matrix for use later:
    sim_file = "$(params["cache"]["cur_dir"])/similarity_matrix.jld2"
    if (isfile(sim_file))
        similarity_matrix = cache_load(sim_file,"similarity_matrix")
    else
        @info "Generating similarity matrix"
        similarity_matrix = NetworkConstruction.coexpression_measure(data_from_dataframe(sample_counts,"data"),params["network_construction"]["coexpression"])
        @info "Saving similarity matrix at $sim_file"
        cache_save(sim_file,"similarity_matrix"=>similarity_matrix)
    end


    ## Adjacency matrix 
    ##update cache
    cache_update("threshold",params["network_construction"]["threshold"])
    cache_update("threshold_method",params["network_construction"]["threshold_method"])
    adj_file = "$(params["cache"]["cur_dir"])/adjacency_matrix.jld2"
    if (isfile(adj_file))
        pre_adj_matrix = cache_load(adj_file,"pre-adj_matrix")
        adj_matrix = cache_load(adj_file,"adjacency_matrix")
    else
        @info "Generating adjacency matrix..."
        if (params["network_construction"]["threshold_method"]=="empirical_dist")
            pre_adj_matrix = NetworkConstruction.empirical_dist_adjacency(similarity_matrix,params["network_construction"]["threshold"])
        elseif (params["network_construction"]["threshold_method"]=="empirical_dist_zero")
            pre_adj_matrix = NetworkConstruction.empirical_dist_zero_adjacency(similarity_matrix,params["network_construction"]["threshold"])
        elseif (params["network_construction"]["threshold_method"]=="hard")
            pre_adj_matrix = NetworkConstruction.adjacency(similarity_matrix,params["network_construction"]["threshold"])
        elseif (params["network_construction"]["threshold_method"]=="top")
            ##TODO setting top x value here for now; should be a parameter, but as an Int rather than Float as threshold param is for other methods
            pre_adj_matrix = NetworkConstruction.top_adjacency(similarity_matrix,10)
        else
            throw(ArgumentError("Threshold method $(params["network_construction"]["threshold_method"]) is not recognised."))

        end
        ##form final adjacency matrix
        adj_matrix = copy(pre_adj_matrix)
        adj_matrix = adj_matrix[:,vec(sum(pre_adj_matrix,dims=1).!=0)]
        adj_matrix = adj_matrix[vec(sum(pre_adj_matrix,dims=2).!=0),:]
        @info "Saving adjacency matrix to $adj_file"
        cache_save(adj_file,["adjacency_matrix"=>adj_matrix,"pre-adj_matrix"=>pre_adj_matrix])
    end 
    #Trim nodes with degree zero
    network_counts = sample_counts[vec(sum(pre_adj_matrix,dims=2).!=0),:]
    #maintain list of vertices in graph
    vertexlist = copy(network_counts[!,:transcript_type])     
    edgelist = NetworkConstruction.edgelist_from_adj(adj_matrix)
    return [adj_matrix, network_counts,vertexlist,edgelist]       
end


function network_visualisation(adj_matrix, network_counts,vertexlist,edgelist)       
    #Synthetic test (just override vertex and edge lists here-- is that ok?)
    if(params["test_name"] == "Synthetic") 
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
end
       
function community_analysis(network_counts,adj_matrix)
    #Network Analysis
    #Type representations 
    ## Community structure
    ##update cache
    anal_dir = "$(params["cache"]["cur_dir"])/communities"
    run(`mkdir -p $(anal_dir)`)
    ##use gene ids here, as they have more chance of getting a GO annotation
    if(params["analysis"]["func_annotate"]==true)
        vertex_gene_names = network_counts[!,:gene_id]
        #community_vertices = GraphletAnalysis.get_community_structure(adj_matrix,vertex_gene_names,"louvain",threejs_plot = true,plot_prefix = "$(params["website"]["website_dir"])/$(params["website"]["page_name"])") 
        community_vertices = GraphletAnalysis.get_community_structure(adj_matrix,vertex_gene_names,"louvain") 
        ## functional annotations of communities
        func_file = "$anal_dir/func_annotations.jld2" 
        if (isfile(func_file))
            functional_annotations = cache_load(func_file,"functional annotations")
        else
            #functional_annotations = GraphletAnalysis.get_functional_annotations(community_vertices,ensembl_version = "75",write_csv = true, csv_dir ="$(params["website"]["website_dir"])/_assets/$(params["website"]["page_name"])/tableinput/")       
            functional_annotations = GraphletAnalysis.get_functional_annotations(community_vertices,ensembl_version = "75")       
            cache_save(func_file,"functional annotations"=>functional_annotations)
        end 
    return [community_vertices,functional_annotations]
    else
        vertex_names = network_counts[!,:transcript_id]
        #community_vertices = GraphletAnalysis.get_community_structure(adj_matrix,vertex_names,"louvain",threejs_plot = true,plot_prefix = "$(params["website"]["website_dir"])/$(params["website"]["page_name"])") 
        community_vertices = GraphletAnalysis.get_community_structure(adj_matrix,vertex_names,"louvain") 
        #find nodes who are not in any community (usually because they are not in connected component).
        community_orphans = findall(x->x==0,in.(1:length(vertexlist),Ref(community_vertices.name)))
        return community_vertices
    end
end

function graphlet_counting(vertexlist,edgelist)

    ##update cache
    anal_dir = "$(params["cache"]["cur_dir"])/graphlets"
    run(`mkdir -p $(anal_dir)`)
    graphlet_file = "$anal_dir/graphlets.jld2" 
    if (isfile(graphlet_file))
        @info "Loading graphlet counts from $anal_dir..."
        graphlet_counts = cache_load(graphlet_file,"graphlets")
        timer = cache_load(graphlet_file,"time")
    else
        @info "Counting graphlets..."
        timer=@elapsed graphlet_counts = GraphletCounting.count_graphlets(vertexlist,edgelist,4,run_method="distributed-old")
        #graphlet_concentrations = concentrate(graphlet_counts) 
        @info "Saving graphlet counts at $anal_dir..."
        ##save the per-edge array as well in case we need it in the future (exp for debugging)
        cache_save(graphlet_file,["graphlets"=>graphlet_counts,"time"=>timer])

    end
    return [graphlet_counts,timer]
end
        
function typed_representations(graphlet_counts,timer,vertexlist,edgelist)
    #method to deduce which run method the null model should use given the run time of graphlets above ($timer)
    if (timer<30)
        null_run = "distributed-short"
    else
        null_run = "distributed-long"
    end

    ##Typed representations
    @info "Looking at typed representations of graphlets..."

    ##update cache
    rep_dir = "$(params["cache"]["cur_dir"])/typed_representations/nullmodel/$(params["analysis"]["null_model_size"])_simulations"
    run(`mkdir -p $(rep_dir)`)
    N=params["analysis"]["null_model_size"]
    rand_graphlets_file = "$rep_dir/rand_graphlets.jld2"
    if (isfile(rand_graphlets_file))
        @info "Loading randomised vertices and graphlet counts from $rep_dir..."
        rand_types_set = cache_load(rand_graphlets_file,"rand vertices")
        rand_graphlet_collection = cache_load(rand_graphlets_file,"rand graphlets")
    else
        ## randomise node types

        #number of randomised graphs
        rand_types_set = [copy(vertexlist) for i in 1:N]
        broadcast(shuffle!,rand_types_set) 
        @info "Counting graphlets on null model" 
        if (null_run=="distributed-short")
            #rand_graphlet_counts = count_graphlets.(rand_types_set,Ref(edgelist),4,run_method="distributed-old")
            rand_graphlet_counts = @showprogress pmap(x->GraphletCounting.count_graphlets(x,edgelist,4,run_method="serial"),rand_types_set,batch_size =10)
        end
        if (null_run=="distributed-long")
            #rand_graphlet_counts = count_graphlets.(rand_types_set,Ref(edgelist),4,run_method="distributed")
            rand_graphlet_counts = @showprogress map(x->count_graphlets(x,edgelist,4,run_method="distributed"),rand_types_set)
        end
        rand_graphlet_collection = vcat(collect.(rand_graphlet_counts)...)
        @info "Saving random graphlet count information at $rep_dir..."
        cache_save(rand_graphlets_file,["rand graphlets"=>rand_graphlet_collection,"rand vertices"=>rand_types_set])
    end


    rand_df = DataFrame(graphlet = broadcast(first,rand_graphlet_collection),value = broadcast(last,rand_graphlet_collection))
    real_df = DataFrame(graphlet = broadcast(first,collect(graphlet_counts)),value = broadcast(last,collect(graphlet_counts)))


    ##convert graphlet_counts dict output to default dictionary, returning 0 for graphlets that don't exist in the real network
    real_dict = DefaultDict(0,graphlet_counts)
    hom_graphlets = unique(last.(split.(unique(real_df[!,:graphlet]),"_")))
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
            het_graphlets = union(first.(split.(real_fil[!,:graphlet],"_4")),first.(split.(rand_fil[!,:graphlet],"_4")))
        elseif (occursin("3",hog))
            het_graphlets = union(first.(split.(real_fil[!,:graphlet],"_3")),first.(split.(real_fil[!,:graphlet],"_3")))
        end
        #store summaries (for TikZ plot)
        summaries = DataFrame()
        for heg in het_graphlets
            rand_fil_fil = filter(:graphlet=>x->x==heg*"_"*hog,rand_df)
            transform!(rand_fil_fil,:value =>ByRow(x-> log(x))=>:log_value)

            ##get summary (for tikZ plot)
            summary = describe(rand_fil_fil[:,3:3],:min,:q25,:median,:q75,:max)
            summary.variable = [heg*"_"*hog]
            append!(summaries,summary)
            ##histogram for each heterogeneous graphlet
            #histogram(rand_fil_fil,:value,:graphlet,"$(params["website"]["website_dir"])/_assets/$(params["website"]["page_name"])/plots/$(heg)_$(hog)_histogram.svg")
            ##log version
            #histogram(rand_fil_fil,:log_value,:graphlet,"$(params["website"]["website_dir"])/_assets/$(params["website"]["page_name"])/plots/$(heg)_$(hog)_log_histogram.svg")
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
        #p = plot(layer(filter(:graphlet=>x->occursin(hog,x),log_real_fil),x = :graphlet,y = :value, Geom.point,color=["count in graph"]),Guide.xticks(label=true),Theme(key_position = :none),Guide.xlabel(nothing),Guide.ylabel("log value"),Guide.yticks(orientation=:vertical),layer(filter(:graphlet=>x->occursin(hog,x),log_rand_fil),x=:graphlet,y=:value,Geom.boxplot(suppress_outliers = true),color=:graphlet));
        #draw(SVG("$(params["website"]["website_dir"])/_assets/$(params["website"]["page_name"])/plots/$(hog)_boxplot.svg",4inch,6inch),p)
        #TeX plot (via PGFPlots) 
        #add real log values to summaries, order from lowest to highest)
        summaries.values = log_real_fil.value
        sort!(summaries,:values)
        #tex_boxplot(summaries[!,Not(:values)],summaries.values,"output/share/$(hog)_boxplot.tex","input",ylabel="")
        ##leave real values attached here, improves tikz layout
        merged_summaries[i] = summaries
        hog_array[i] = hog_df
        hog_array_under[i] = hog_df_under
    end

    #merged boxplots
    #tex_merged_boxplot(merged_summaries,"output/share/merged_boxplot.tex","input",ylabel = "log value")

    ## find significant graphlets
    sig_graphlets = vcat(filter.(:p_value=>x->x<0.05,hog_array)...)
    insig_graphlets = vcat(filter.(:p_value=>x->x<0.05,hog_array_under)...)

    #save in output cache
    NetworkConstruction.html_table_maker(sig_graphlets,"$rep_dir/sig_type_representations.html",imgs=sig_graphlets.Graphlet)                          
    NetworkConstruction.html_table_maker(insig_graphlets,"$rep_dir/insig_type_representations.html",imgs=sig_graphlets.Graphlet)                              
    #save for website version
    #NetworkConstruction.html_table_maker(sig_graphlets,"$(params["website"]["website_dir"])/_assets/$(params["website"]["page_name"])/sig_type_representations.html",imgs=sig_graphlets.Graphlet,figpath="../figs/")
    #NetworkConstruction.html_table_maker(insig_graphlets,"$(params["website"]["website_dir"])/_assets/$(params["website"]["page_name"])/insig_type_representations.html",imgs=insig_graphlets.Graphlet,figpath="../figs/")  
    ##look at edge types in randomised networks
    real_type_edgecounts = countmap(splat(tuple).(sort.(eachrow(hcat(map(x->vertexlist[x],first.(edgelist)),map(x->vertexlist[x],last.(edgelist)))))))
    rand_types_edgecounts = map(y->(countmap(splat(tuple).(sort.(eachrow(hcat(map(x->y[x],first.(edgelist)),map(x->y[x],last.(edgelist)))))))),rand_types_set)
    rand_edge_collection = vcat(collect.(rand_types_edgecounts)...)
    rand_edge_df = DataFrame(graphlet = broadcast(first,rand_edge_collection),value = broadcast(last,rand_edge_collection))
    random_edges = DataFrame()
    for t in unique(rand_edge_df[!,:graphlet])
        rand_vals = filter(:graphlet=>x->x==t,rand_edge_df)[!,:value]
        rand_exp = sum(rand_vals)/N
        real_obs = real_type_edgecounts[t]
        append!(random_edges,DataFrame(Graphlet = first(t)*"_"*last(t)*"_edge", Expected = rand_exp,Observed = real_obs))       
    end

    #pretty_table(random_edges,backend=:html,standalone = false)
end 
                #@time motif_counts = find_motifs(edgelist,"hetero_rewire",100, typed = true, typelist = vec(vertexlist),plotfile="$cache_dir/motif_detection.svg",graphlet_size = 4)



#               ### Validation steps
                #val_dir = "$anal_dir/validation"
                #run(`mkdir -p $(val_dir)`)
#               # looking at identified significant graphlets and seeing if they check out biologically
#               # the most taxing step is to identify the graphlets that are coincident in some way to KEGG pathways. We cache these coincidents as a dataframe (using CSV instead of JLD)  
                #
function coincident_analysis(adj_matrix,network_counts,vertexlist,edgelist)
    vertex_names = network_counts[!,:transcript_id]
    #Coincident analysis
    coinc_dir = "$(params["cache"]["cur_dir"])/coincidents"
    run(`mkdir -p $(coinc_dir)`)
    #get baseline entrez and kegg info about transcripts
    kegg_file = "$(coinc_dir)/kegg_info.jld2"
    if (isfile(kegg_file))
        @info "Loading KEGG info from $kegg_file..."
        entrez_id_vector = cache_load(kegg_file,"entrez_id_vector")
        candidates = cache_load(kegg_file,"candidates")
        top_terms = cache_load(kegg_file,"top_terms")
    else

        @info "getting KEGG info..."

        entrez_id_vector, candidates,top_terms = GraphletAnalysis.get_KEGG_pathways(vertex_names,"transcripts")
        cache_save(kegg_file,["entrez_id_vector"=>entrez_id_vector, "candidates"=>candidates,"top_terms"=>top_terms ])
    end

    candidate_pathways = collect(keys(candidates))

    coincidents_file ="$coinc_dir/coincidents.csv"
    if (isfile(coincidents_file))
        @info "Loading coincidents dataframe from $coinc_dir..."
        Coincidents = CSV.read(coincidents_file,DataFrame)
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
        Coincidents = GraphletAnalysis.graphlet_coincidences(vertexlist,vertex_names,"transcripts",adj_matrix,entrez_id_vector,candidates)
        @info "Saving coincidents at $coinc_dir..."
        CSV.write(coincidents_file,Coincidents)

    end


    #find which types are excluded in general, and then only as a cause of having no Entrez id
    Coincidents.excluded = [Coincidents.Transcript_type[i][Coincidents.Inclusion[i].==0] for i in 1:size(Coincidents)[1]]
    Coincidents.excluded_Entrez = [Coincidents.Transcript_type[i][Coincidents.Entrez[i].==0] for i in 1:size(Coincidents)[1]]
    #find only those coincidents that involve non-coding transcripts
    Coincidents_noncoding = Coincidents[findall(x-> "noncoding" in x, Coincidents.Transcript_type),:]

    #Nonuniformity test: finds graphlets where the included nodes differ across different pathways (sampling for now for speed)
    graphlet = "4-star"
    nonuniforms = []
    for y in filter(p->last(p)>1,countmap(filter(:Hom_graphlet=> x-> x == graphlet, Coincidents[1:199000,:]).Vertices))
        test = sum(filter(:Vertices => x-> x == first(y),filter(:Hom_graphlet => x-> x == graphlet,Coincidents))[!,8])/last(y)
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
    #orbitsigs_file ="$val_dir/orbitsigs.jld"
    #if (isfile(orbitsigs_file))
    #        @info "Loading orbitsigs from $val_dir..."

    #collect dataframe for each node in this array
    #per_node_significance = Array{DataFrame,1}(undef,length(vertexlist))
    #choose just one order of graphlets (3 or 4)
    sub_Coincidents = filter(:Hom_graphlet=>x->occursin("4-",x),Coincidents)
    orbit_sigs_file = "$coinc_dir/orbit_sigs.jld2" 
    ## select whether we are looking at "detailed" or "collated" significance
    orbit_sigs_method ="detailed"
    
    if (isfile(orbit_sigs_file))
        @info "Loading orbit significance dataframe from $coinc_dir..."
        orbit_sigs = cache_load(orbit_sigs_file,"orbit_sigs")
    else

        @info "Getting orbit significance statistics..."
        #table to showing whether each node (row) is included in each pathway (column)
        inkey = hcat([ in.(1:length(vertexlist),Ref(candidates[p])) for p in candidate_pathways ]...)
        if(orbit_sigs_method == "collated")
            orbit_sigs = @showprogress map(x->GraphletAnalysis.pernode_significance(x,sub_Coincidents,candidate_pathways,inkey[x,:]),1:length(vertexlist))
        end
        if(orbit_sigs_method == "detailed")
            orbit_sigs = @showprogress map(x->GraphletAnalysis.pernode_significance_detail(x,sub_Coincidents,candidate_pathways,inkey[x,:]),1:length(vertexlist))
        end
        cache_save(orbit_sigs_file,"orbit_sigs"=>orbit_sigs)
    end

    ## now compare the significance profile of those nodes that are not attached to a pathway to the average pathway profile of known pathway nodes
    ##convert to array form for comparisons
    orbit_sigs_array = map(x->Array(x[!,2:end]),orbit_sigs)
    # this table stores the average for each pathway (of nodes that are known to be in the pathway)
    significance_bars = zeros(length(keys(candidates)),size(orbit_sigs_array[1])[2])
    for (i,c) in enumerate(keys(candidates))
        significance_bars[i,:] = (sum(map(x->Array(x[!,2:end]),orbit_sigs[candidates[c]]))./length(vertexlist))[i,:]
    end
    ##now compare bars against profiles of non-pathway nodes
    ## choose a subset of nodes to look at. Can be boolean BitArray (with length equal to all nodes) or a specific list of nodes 
    subset = entrez_id_vector.==0
    putative_pathways = Array{Array{String,1}}(undef,length(orbit_sigs_array[subset]))
    for (i,t) in enumerate(orbit_sigs_array[subset])
        putative_pathways[i] = collect(keys(candidates))[vec((sum(t.>significance_bars,dims=2).>2))]
    end

    # Gadfly beeswarm visualisation:
    output_dir = "$cwd/output/plots/orbit_significance_$(orbit_sigs_method)/"
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
    #TODO remove this and just use low_filter method below to more organically acheive same thing
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
    zero_orbit_sigs_array = map(x->Array(x[!,2:end]),zero_orbit_sigs)
    zero_candidates = Dict(Pair.(zero_candidate_pathways,[candidates[x] for x in zero_candidate_pathways]))

    #uniqueness of pathway contributors (concerns of too much overlap) 
    sig_pathway_occurences = countmap(vcat([candidates[x] for x in zero_candidate_pathways]...))
    m = max(collect(values(sig_pathway_occurences))...) 
    #m = 8
    supersharers = first.(filter(x->last(x)==m,collect(sig_pathway_occurences)))
    #for these supersharers, find the set of pathways they are involved in
    supersharer_pathways = GraphletAnalysis.pathways_per_node_dict(supersharers,zero_candidates)
    in_group = collect(keys(countmap(vcat(collect(values(supersharer_pathways))...))))
    not_in_group = zero_candidate_pathways[.!(in.(zero_candidate_pathways,Ref(collect(keys(countmap(vcat(collect(values(supersharer_pathways))...)))))))]
    countmap(collect(values(supersharer_pathways)))
    #do we need to rule out pathways dominated by supersharers? TODO


    #ecdfs
    #store ecdf functions in table
    ecdf_table = Array{ECDF,2}(undef,length(zero_candidate_pathways),last_col)
    for (i,p) in enumerate(zero_candidate_pathways)
        for j in 1:last_col
            ecdf_table[i,j] = ecdf(map(x->x[i,j],zero_orbit_sigs_array))
        end
    end
    #first check candidate probabilities across all categories
    known_pathway_probs = Array{Array{Float64,2},1}(undef,length(zero_candidate_pathways))
    for (i,c) in enumerate(zero_candidate_pathways)
        known_pathway_probs[i] = hcat([map(ecdf_table[i,j],map(x->x[i,j],orbit_sigs_array[zero_candidates[c]])) for j in 1:last_col]...)
    end

    #define orbit names here for now
    orbit_names = names(orbit_sigs[1])[2:last_col+1]

    #collect known pathway vectors for corresponding known pathway nodes
    known_pathway_dfs = Array{DataFrame,1}(undef,length(zero_candidate_pathways))
    for (i,p) in enumerate(zero_candidate_pathways)
        subset = zero_candidates[p]
        df_build = DataFrame(hcat(map(x->x[i,:],zero_orbit_sigs_array[subset])...)',orbit_names)
        insertcols!(df_build,1,:shared=>[sig_pathway_occurences[x] for x in subset].-1)
        insertcols!(df_build,1,:transcript_id=>subset)
        known_pathway_dfs[i] = df_build 
        #print("Known nodes in pathway $p are in this many other pathways:\n")
        #print([sig_pathway_occurences[x] for x in subset].-1)
        #print("\n")
    end
    known_pathway_arrays = map(x->Array(x[!,3:end]),known_pathway_dfs)
    #known ecdfs
    #calculate ecdfs just for known pathway node values (each row corresponds to a pathway, each column to an orbit) 
    known_ecdf_table = Array{ECDF,2}(undef,length(zero_candidate_pathways),last_col)
    for (i,p) in enumerate(zero_candidate_pathways)
        known_ecdf_table[i,:] = ecdf.(eachcol(known_pathway_arrays[i]))
    end

    ##compare each unknown node to known ecdfs
    #set threshold for which to check on i.e. is node orbit count higher than the top ~thresh~% of known nodes
    thresh = 0.05
    ub = 1 -thresh
    #first check if there are any orbits for any pathways that have such low representation that 0 count exceed threshold. we will filter these out of analysis
    low_filter = .!(map(x->x(0),known_ecdf_table).>ub)
    #for each node, map each known ecdf to the corresponding orbit count and check against threshold, factoring in the low filter
    unknown_ecdf_comparison = map(y->(reshape(map(x->known_ecdf_table[x](y[x]),1:length(known_ecdf_table)),size(known_ecdf_table)[1],size(known_ecdf_table)[2]).>ub).*low_filter,zero_orbit_sigs_array)
    #determine a node as significantly linked to a pathway if at least half its orbit counts are significant
    sig_check = map(x->(sum(x,dims=2).>(last_col/2)),unknown_ecdf_comparison)
    sig_nodes= findall(x->sum(x)>0,sig_check)
    sig_nodes_dict = Dict(Pair.(sig_nodes,map(x->zero_candidate_pathways[vec(x)],sig_check[sig_nodes])))

    #shape plots into a grid
    ncols = 6
    dims = fldmod(length(zero_candidate_pathways),ncols)
    plots = Array{Union{Plot,Context},2}(undef,dims[1]+(dims[2]>0),ncols)
    #plots = Array{Union{Plot,Context},1}(undef,length(keys(candidates)))
    for (i,p) in enumerate(zero_candidate_pathways)
        #For total ecdfs:
        #find max over all measures for pathway
        #                    m = max([map(x->x[i,1],orbit_sigs_array)...,map(x->x[i,2],orbit_sigs_array)..., map(x->x[i,3],orbit_sigs_array)...]...)
        #                    plots[i] = plot([layer(x->ecdf(map(x->x[i,j],orbit_sigs_array))(x),0,m,color=[j]) for j in 1:last_col]...,
        #                              Scale.color_discrete_manual("orange", "green", "purple"),
        #                              Guide.title(p),
        #                              Guide.xlabel("count"),
        #                              Theme(major_label_font_size=4pt,key_position=:none));
        #                              #Guide.colorkey(title="orbit position"),
        #                              #Guide.title(p));
        #                              #
        #for known ecdfs only:
        m = max(known_pathway_arrays[i]...)
        plots[i] = plot([layer(x->known_ecdf_table[i,j](x),0,m,color=[orbit_names[j]]) for j in 1:last_col]...,
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
    draw(SVG("$(output_dir)/known_ecdfs.svg",30cm,20cm),gridstack(plots))


    ## make a table to analyse different candidate levels for coincident graphlets


    ##function to get the n permutations of a set xs
    all_perm(xs, n) = vec(map(collect, Iterators.product(ntuple(_ -> xs, n)...)))
    ## method to find all orbit permutations of a type set
    set = ["coding","noncoding"]

    # specify orbit classes
    orbit = [1,1,1,1]
    order = length(orbit)
    ##get all possible combinations from set (before sorting)
    base = collect.(vcat(collect(Iterators.product(repeat([set],order)...))...))
    #sort via each orbit equivalence class 
    for o in unique(orbit)
        template = orbit.==o
        for i in base
            i[template] = sort(i[template])
        end
    end
    #find unique permutations
    gs = unique(base)  

    return sig_nodes_dict
end                    



