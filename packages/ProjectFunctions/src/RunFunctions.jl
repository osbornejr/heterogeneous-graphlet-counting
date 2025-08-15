#using Distributed, JLD2, CSV
#using StatsBase,Gadfly,Compose,DataFrames,YAML,Dates, Graphs, Colors, Random, Distributions, ProgressMeter
using Pkg,ProgressMeter,DataFrames,YAML,Distributed,JLD2,CSV,StatsBase,Random,Graphs,Dates,Colors,Gadfly,Compose,DataStructures,CategoricalArrays,GraphPlot



function run_all(config_file::String)
    ##TODO move the adding/removing of workers to align with the individual functions that require distributed framework (adds sometimes unnessesary load time to start of run otherwise). Turning off the global initialisation here for now   
#    if (length(workers())!=Threads.nthreads())
#        #Pkg.resolve()
#       #Pkg.precompile()
#        @info "Setting up worker processes"
#        distributed_setup(:ProjectFunctions,:GraphletCounting,:GraphletAnalysis,:NetworkConstruction)
#    end

    @info "Loading parameters"
    load_config(config_file)
    @info "Loading raw counts"
    raw_counts = get_input_data()
    @info "Preprocessing raw counts"
    processed_counts = data_preprocessing(raw_counts)
    if (params["network_construction"]["construct"] == true)
        @info "Constructing network"
        adj_matrix,network_counts,vertexlist,edgelist = network_construction(processed_counts)
        anal_flag = 0
    else
        #do not analyse if no network constructed
        anal_flag = 1
    end

    ## For analysis step, we need to run for both real and synthetic networks
    #this means we will need to reload the cache dirs for each run
    #run first on synthetic network, as that is what will be loaded in if "synthetic" parameter is true
    # we then reload network_construction with "synthetic" updated to false
    #
    while (anal_flag == 0)  
        
        #only do bio val on real network 
        if (params["analysis"]["biological_validation"] * !params["network_construction"]["synthetic"] == true  )
            @info "Conducting baseline network validation"
            kegg_terms,go_terms = biological_validation(network_counts)
        end


        if (params["analysis"]["community_detection"] == true)
            @info "Finding communities"
            com_anal = community_analysis(network_counts,adj_matrix)
        end
    
        if (params["analysis"]["wgcna"] == true)
            @info "Finding WGCNA comparison network"
            wgcna_anal = wgcna_analysis(processed_counts)
        end

        if (params["analysis"]["graphlet_counting"] == true)
            @info "Counting graphlets"
            graphlet_counts,timer = graphlet_counting(vertexlist,edgelist)
        end

        if (params["analysis"]["typed_representation"] == true)
            @info "Comparing typed graphlet representations"
            typed_anal = typed_representations(graphlet_counts,timer,vertexlist,edgelist)
        end
        
        #note that coincident analysis only makes sense on real data TODO split out enumeration and coincident analysis to allow synthetic enumeration
        if (params["analysis"]["coincident_graphlets"] * !params["network_construction"]["synthetic"] == true  )
            @info "Conducting coincident graphlet analysis"
            coinc_graphlets = coincident_graphlets(network_counts,vertexlist,edgelist)
        end
        
        ##reload analysis if necessary
        if (params["network_construction"]["synthetic"] == true)
            #synthetic run complete, update cache directories for real network
            params["network_construction"]["synthetic"] = false
            cache_setup()
            @info "Synthetic network analysis complete, loading in real network for analysis"
            adj_matrix,network_counts,vertexlist,edgelist = network_construction(processed_counts)
        else
            anal_flag = 1
        end
    end


end 

function load_config(config_file::String)
        #setup parameters to be read.writable for all functions in this module
        global params = YAML.load_file(config_file)
        cache_setup()
end
export load_config
export params

#exposed function to allow just preprocessing run on given configs settings
function data_preprocessing(config_file::String)
    @info "Loading parameters"
    load_config(config_file)
    @info "Loading raw counts"
    raw_counts = get_input_data()
    @info "Preprocessing raw counts"
    processed_counts = data_preprocessing(raw_counts)
end

function data_preprocessing(raw_counts::DataFrame)
    
    #TODO may need to add some kind of deduplication step!
    
    ##Round: round reads to significant figures (TODO only applicable for FPKM/expected counts? Maybe need a switch in config for this)

    file = "$(params["cache"]["round_dir"])/round_counts.jld2"

        
    if(isfile(file))
        round_counts = cache_load(file,"round counts")
    else

        round_counts = DataPreprocessing.round_raw_counts(raw_counts,params["data_preprocessing"]["roundsig"])
        @info "Saving round counts to $file"
        cache_save(file,"round counts"=>round_counts)
    end

    ### VST- log transform data before applying any other processing.

    file = "$(params["cache"]["vst_dir"])/vst_counts.jld2"

        
    if(isfile(file))
        vst_counts = cache_load(file,"vst counts")
    else
        if (params["data_preprocessing"]["vst"] == true) 
            vst_counts = DataPreprocessing.log_counts(round_counts)
        else
            vst_counts = round_counts
        end

        @info "Saving vst counts to $file"
        cache_save(file,"vst counts"=>vst_counts)
    end

    ## Clean - remove transcripts with total counts across all samples less than Cut
    file = "$(params["cache"]["clean_dir"])/clean_counts.jld2"

    
    if(isfile(file))
        clean_counts = cache_load(file,"clean counts")
    else

        clean_counts = DataPreprocessing.clean_round_counts(vst_counts,params["data_preprocessing"]["expression_cutoff"],params["data_preprocessing"]["minreq"],method=params["data_preprocessing"]["clean_method"])
        @info "Saving clean counts to $file"
        cache_save(file,"clean counts"=>clean_counts)
    end


    ### Normalisation
    ##update cache
    file = "$(params["cache"]["norm_dir"])/norm_counts.jld2"

    if(isfile(file))
        norm_counts = cache_load(file,"norm counts")
    else
        @info "Saving norm counts to $file"
        norm_counts = DataPreprocessing.normalise_clean_counts(clean_counts,params["data_preprocessing"]["norm_method"])
        cache_save(file,"norm counts"=>norm_counts)
    end

    ### High contrast pruning
    ##update cache
    file = "$(params["cache"]["prune_dir"])/prune_counts.jld2"
    if(isfile(file))
        prune_counts = cache_load(file,"prune counts")
    else
        if (params["data_preprocessing"]["prune"] == true) 
           #TODO needs to work for human data with mroe than two options for a condition 
           prune_counts = DataPreprocessing.prune_normalised_counts(norm_counts,params["conditions"][params["data_preprocessing"]["prune_condition"]]["samples_1"],params["conditions"][params["data_preprocessing"]["prune_condition"]]["samples_2"],params["data_preprocessing"]["prune_cutoff"],params["data_preprocessing"]["prune_strictness"])
        else
        prune_counts = norm_counts
        end

        @info "Saving pruned counts to $file"
        cache_save(file,"prune counts"=>prune_counts)
    end


    ##Sampling for most variable transcripts
    ##update cache
    file = "$(params["cache"]["sampling_dir"])/sample_counts.jld2"

    if(isfile(file))
        sample_counts = cache_load(file,"sample counts")
    else
        @info "Saving sample counts to $file"
        sample_counts = DataPreprocessing.sample_norm_counts(prune_counts,params["data_preprocessing"]["variance_percent"],method=params["data_preprocessing"]["sample_method"])
        cache_save(file,"sample counts"=>sample_counts)
    end

    return sample_counts
end 
export data_preprocessing


function network_construction(sample_counts::DataFrame)

    ##Network construction
    ##Measure of coexpression
    ##update cache


    #similarity_matrix=mutual_information(data)
    ## file to cache similarity matrix for use later:
    sim_file = "$(params["cache"]["similarity_dir"])/similarity_matrix.jld2"
     if (isfile(sim_file))
        similarity_matrix = cache_load(sim_file,"similarity_matrix")
    else
        @info "Generating similarity matrix"
        #Setup workers
        distributed_setup(:ProjectFunctions,:GraphletCounting,:GraphletAnalysis,:NetworkConstruction)
        
        similarity_matrix = coexpression_measure(data_from_dataframe(sample_counts,"data"),params["network_construction"]["coexpression"])
        @info "Saving similarity matrix at $sim_file"
        cache_save(sim_file,"similarity_matrix"=>similarity_matrix)
        ##cleanup distributed
        rmprocs(workers())
        cleaner()
    end


    ## Adjacency matrix 
    ##update cache
    adj_file = "$(params["cache"]["adjacency_dir"])/adjacency_matrix.jld2"
    if (isfile(adj_file))
        pre_adj_matrix = cache_load(adj_file,"pre-adj_matrix")
        adj_matrix = cache_load(adj_file,"adjacency_matrix")
        components =cache_load(adj_file,"components")

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
        ##form final adjacency matrix (remove zero nodes)
        adj_matrix = copy(pre_adj_matrix)
        
        ##TODO HACK: remove dodgy second component from network (kinda assumes there are no zero nodes)
        #get components
        components = NetworkConstruction.network_components(pre_adj_matrix)
        #get largest component
        largest = findmax(length.(components))[2]
        adj_matrix = pre_adj_matrix[components[largest],components[largest]]
        
        #old method-- just get all nonzero components
        #adj_matrix = adj_matrix[:,vec(sum(pre_adj_matrix,dims=1).!=0)]
        #adj_matrix = adj_matrix[vec(sum(pre_adj_matrix,dims=1).!=0),:]



        @info "Saving adjacency matrix to $adj_file"
        cache_save(adj_file,["adjacency_matrix"=>adj_matrix,"pre-adj_matrix"=>pre_adj_matrix,"components"=>components])
    end 


    

    #get largest component
    largest = findmax(length.(components))[2]

    network_counts = sample_counts[components[largest],:]
    #or to just get all non-zero components
    #Trim nodes with degree zero (and TODO HACK remove component 2 as well!
    #network_counts = sample_counts[vec(sum(pre_adj_matrix,dims=2).!=0),:]
    #maintain list of vertices in graph
    vertexlist = copy(network_counts[!,:transcript_type])     
    edgelist = NetworkConstruction.edgelist_from_adj(adj_matrix)

    if(params["network_construction"]["synthetic"] == true)
        ##slightly hacky way to do this still (TODO does the vertexlist (types) need to be randomised? it is atm)
        @info "Switching to synthetic network"
        edgelist,vertexlist = NetworkConstruction.synthetic_network(vertexlist,edgelist)
        adj_matrix = NetworkConstruction.adj_from_edgelist(edgelist)
    end
    return [adj_matrix, network_counts,vertexlist,edgelist]       
end
export network_construction

##catch all function to call desired coexpression_measure 
function coexpression_measure(data::Union{AbstractDataFrame,AbstractArray},method::String)
    ## wrapper for all possible different coexpression measure methods
    ## convert to array if necessary
    if(typeof(data)<:AbstractDataFrame)
        data = DataPreprocessing.data_from_dataframe(data,"data")
    end 

    if (method =="pearson")
        return cor(data')
    end

    if (method=="spearman")
        return corspearman(data')
    end

    if (method== "kendall")
        return corkendall(data')
    end

    if (method== "ridge_partialcor")
        return NetworkConstruction.ridge_partialcor(data',ridge=1e-3)
    end

    if (method=="mutual_information")
        
        ##determine nbins from config
        nvars,nvals = size(data)
        if params["network_construction"]["nbins"] == "sqrt_n"
            nbs = Int(round(sqrt(nvals)))
        elseif params["network_construction"]["nbins"] == "n_over_2"
            nbs = Int(round(nvals/2))
        end
        return NetworkConstruction.mutual_information(data; discretizer = "uniform_width", estimator = "maximum_likelihood", mi_base = 2,nbins = nbs )
    end
    if (method=="pidc")
        ##determine nbins from config
        nvars,nvals = size(data)
        if params["network_construction"]["nbins"] == "sqrt_n"
            nbs = Int(round(sqrt(nvals)))
        elseif params["network_construction"]["nbins"] == "n_over_2"
            nbs = Int(round(nvals/2))
        end
        return NetworkConstruction.partial_information_decomposition(data; discretizer = "uniform_width", estimator = "maximum_likelihood", mi_base = 2,distributed = true, nbins = nbs)
    end
    if (method=="pcit")
        return NetworkConstruction.pcit(data)
    end
end

function network_visualisation(adj_matrix, network_counts,vertexlist,edgelist)       

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
    ## Community structure
    ##update cache
    anal_dir = params["cache"]["community_dir"]    
    run(`mkdir -p $(anal_dir)`)

    ##use gene ids here, as they have more chance of getting a GO annotation
    if(params["analysis"]["func_annotate"]==true)
        vertex_gene_names = network_counts[!,:gene_id]
        # if (params["website"]["website"] == true)
        #     community_vertices = GraphletAnalysis.get_community_structure(adj_matrix,vertex_gene_names,"louvain",threejs_plot = true,plot_prefix = "$(params["website"]["website_dir"])/$(params["website"]["page_name"])") 
        # else
        communities_file = "$anal_dir/communities.jld2" 
        if (isfile(communities_file))
            community_vertices = cache_load(communities_file,"communities")
        else
            community_vertices = GraphletAnalysis.get_community_structure(adj_matrix,vertex_gene_names,"louvain") 
            #save communities
            cache_save(communities_file,"communities"=>community_vertices)
        end
        ## functional annotations of communities
        func_file = "$anal_dir/func_annotations.jld2" 
        if (isfile(func_file))
            functional_annotations = cache_load(func_file,"functional annotations")
        else
            #functional_annotations = GraphletAnalysis.get_functional_annotations(community_vertices,ensembl_version = "75",write_csv = true, csv_dir ="$(params["website"]["website_dir"])/_assets/$(params["website"]["page_name"])/tableinput/")       
            #functional_annotations = GraphletAnalysis.get_functional_annotations(community_vertices,ensembl_version = "75")       
            functional_annotations = GraphletAnalysis.get_functional_annotations(community_vertices)       
            cache_save(func_file,"functional annotations"=>functional_annotations)
        end 
        return [community_vertices,functional_annotations]

    else
        vertex_names = network_counts[!,:transcript_id]
        # if (params["website"]["website"] == true)
        #     community_vertices = GraphletAnalysis.get_community_structure(adj_matrix,vertex_names,"louvain",threejs_plot = true,plot_prefix = "$(params["website"]["website_dir"])/$(params["website"]["page_name"])") 
        # else
        communities_file = "$anal_dir/communities.jld2" 
        if (isfile(communities_file))
            community_vertices = cache_load(communities_file,"communities")
        else
            community_vertices = GraphletAnalysis.get_community_structure(adj_matrix,vertex_names,"louvain") 
            #save communities
            

            cache_save(communities_file,"communities"=>community_vertices)
        end
        return community_vertices
    end
    #find nodes who are not in any community (usually because they are not in connected component).
    #community_orphans = findall(x->x==0,in.(1:length(vertexlist),Ref(community_vertices.name)))

end
export community_analysis

function wgcna_analysis(processed_counts)
    ##WGCNA needs to set up its own network (on the same processed counts) so we output both the network and the derived community structure here
    #TODO integrate into community step somehow? subset or separate?
    
    wgcna_dir = params["cache"]["wgcna_dir"]    
    run(`mkdir -p $(wgcna_dir)`)
    wgcna_file = "$wgcna_dir/wgcna.jld2"
    if (isfile(wgcna_file))
        @info "Loading WGCNA information from $wgcna_file..."
        wgcna_network = cache_load(wgcna_file,"adjacency")
        wgcna_comm = cache_load(wgcna_file,"communities")
    else
        @info "Running WGCNA..." 
        wgcna_network,groups,colours = NetworkConstruction.wgcna(data_from_dataframe(processed_counts),processed_counts.transcript_type)
        wgcna_comm = DataFrame(label = processed_counts.transcript_id,
                               group = groups,
                               color = "#".*(colours)
                              )

        cache_save(wgcna_file,["adjacency"=>wgcna_network,"communities"=>wgcna_comm])
    end
    return [wgcna_network,wgcna_comm]
end
export wgcna_analysis

function graphlet_counting(vertexlist,edgelist)

    ##update cache
    anal_dir = params["cache"]["graphlet_counting_dir"]    
    run(`mkdir -p $(anal_dir)`)
    graphlet_file = "$anal_dir/graphlet-counting.jld2" 
    if (isfile(graphlet_file))
        @info "Loading graphlet counts from $anal_dir..."
        graphlet_counts = cache_load(graphlet_file,"graphlets")
        timer = cache_load(graphlet_file,"time")
    else
        @info "Counting graphlets..."
        #Setup workers
        distributed_setup(:ProjectFunctions,:GraphletCounting,:GraphletAnalysis,:NetworkConstruction)

        timer=@elapsed graphlet_counts = GraphletCounting.count_graphlets(vertexlist,edgelist,params["analysis"]["graphlet_size"],run_method="distributed-old")
        #graphlet_concentrations = concentrate(graphlet_counts) 
        @info "Saving graphlet counts at $anal_dir..."
        ##save the per-edge array as well in case we need it in the future (exp for debugging)
        cache_save(graphlet_file,["graphlets"=>graphlet_counts,"time"=>timer])
        ##cleanup distributed
        rmprocs(workers())
        cleaner()

    end
    return [graphlet_counts,timer]
end 
export graphlet_counting

function typed_representations(graphlet_counts,timer,vertexlist,edgelist)
    #TODO rewrite this function so that the guts of it is in appropriate package
    #method to deduce which run method the null model should use given the run time of graphlets above ($timer)
    if (timer<30)
        null_run = "distributed-short"
    else
        null_run = "distributed-long"
    end

    ##Typed representations
    @info "Looking at typed representations of graphlets..."

    ##update cache
    rep_dir = params["cache"]["rep_dir"]    
    run(`mkdir -p $(rep_dir)`)
    N=params["analysis"]["null_model_size"]
    rand_graphlets_file = "$rep_dir/rand_graphlets.jld2"
    if (isfile(rand_graphlets_file))
        @info "Loading randomised vertices and graphlet counts from $rep_dir..."
        rand_types_set = cache_load(rand_graphlets_file,"rand vertices")
        rand_graphlet_collection = cache_load(rand_graphlets_file,"rand graphlets")
    else
        ## randomise node types

        #Setup workers
        distributed_setup(:ProjectFunctions,:GraphletCounting,:GraphletAnalysis,:NetworkConstruction,nthreads=24)
        
        #number of randomised graphs
        rand_types_set = [copy(vertexlist) for i in 1:N]
        broadcast(shuffle!,rand_types_set) 
        @info "Counting graphlets on null model" 
        order = params["analysis"]["graphlet_size"]
        if (null_run=="distributed-short")
            #rand_graphlet_counts = count_graphlets.(rand_types_set,Ref(edgelist),4,run_method="distributed-old")
            rand_graphlet_counts = @showprogress pmap(x->GraphletCounting.count_graphlets(x,edgelist,order,run_method="serial"),rand_types_set,batch_size =10)
        end
        if (null_run=="distributed-long")
            #rand_graphlet_counts = count_graphlets.(rand_types_set,Ref(edgelist),4,run_method="distributed")
            rand_graphlet_counts = @showprogress map(x->GraphletCounting.count_graphlets(x,edgelist,order,run_method="distributed"),rand_types_set)
        end
        rand_graphlet_collection = vcat(collect.(rand_graphlet_counts)...)
        @info "Saving random graphlet count information at $rep_dir..."
        cache_save(rand_graphlets_file,["rand graphlets"=>rand_graphlet_collection,"rand vertices"=>rand_types_set])
        ##cleanup distributed
        rmprocs(workers())
        cleaner()
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
    ##get all null model network values individually
    #null_model_arrays = Array{Array}(undef,length(hom_graphlets))

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
        ##if a het graphlet does not occur in real network, add it with a 0 count
        for l in setdiff(het_graphlets,replace.(real_fil.graphlet,"_$hog"=>""))
            append!(real_fil,DataFrame(graphlet=l,value=0))  
        end
        

        #store summaries (for TikZ plot)
        summaries = DataFrame()
        
        #null_model_array = Matrix{N,length(het_graphlets)}

        for heg in het_graphlets
            rand_fil_fil = filter(:graphlet=>x->x==heg*"_"*hog,rand_df)
            transform!(rand_fil_fil,:value =>ByRow(x-> log(x))=>:log_value)
            
            ##get summary (for tikZ plot)
            summary = describe(rand_fil_fil[:,3:3],:min,:q25,:median,:q75,:max)
            summary.variable = [heg*"_"*hog]
            append!(summaries,summary)
            ##histogram for each heterogeneous graphlet
            #
            if (params["website"]["website"] == true)
                DataPreprocessing.histogram(rand_fil_fil,:value,:graphlet,"$(params["website"]["website_dir"])/_assets/$(params["website"]["page_name"])/plots/$(heg)_$(hog)_histogram.svg")
                ##log version
                #histogram(rand_fil_fil,:log_value,:graphlet,"$(params["website"]["website_dir"])/_assets/$(params["website"]["page_name"])/plots/$(heg)_$(hog)_log_histogram.svg")
            end
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
        if (params["website"]["website"] == true)
            #SVG plot (for web)
            p = plot(layer(filter(:graphlet=>x->occursin(hog,x),log_real_fil),x = :graphlet,y = :value, Geom.point,color=["count in graph"]),Guide.xticks(label=true),Theme(key_position = :none),Guide.xlabel(nothing),Guide.ylabel("log value"),Guide.yticks(orientation=:vertical),layer(filter(:graphlet=>x->occursin(hog,x),log_rand_fil),x=:graphlet,y=:value,Geom.boxplot(suppress_outliers = true),color=:graphlet));
            draw(SVG("$(params["website"]["website_dir"])/_assets/$(params["website"]["page_name"])/plots/$(hog)_boxplot.svg",4inch,6inch),p)
        end
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
    #NetworkConstruction.html_table_maker(sig_graphlets,"$rep_dir/sig_type_representations.html",imgs=sig_graphlets.Graphlet,figpath = "$cwd/website/figs")                          
    #NetworkConstruction.html_table_maker(insig_graphlets,"$rep_dir/insig_type_representations.html",imgs=sig_graphlets.Graphlet,figpath = "$cwd/website/figs")                          
    #
    ##For paper 
    #TODO(shift these to proper place (Visualisation?)
    #Over and underpresented graphlets table
    #NetworkConstruction.tex_table_maker(sig_graphlets,"output/share/overrepresented_graphlets.tex")
    #NetworkConstruction.tex_table_maker(insig_graphlets,"output/share/underrepresented_graphlets.tex")
    ##boxplots
    

 #   for s in merged_summaries
 #       ##get homogeneous graphlet from summary
 #       hog = split(s.variable[1],"_")[end]
 #       NetworkConstruction.tex_boxplot(s[!,Not(:values)],s.values,"output/share/$(hog)_boxplot.tex","input",ylabel="")
 #   end
 #   #merged boxplots
 #   NetworkConstruction.tex_merged_boxplot(merged_summaries,"output/share/merged_boxplot.tex","input",ylabel = "log value")


    if (params["website"]["website"] == true)
        #save for website version
        NetworkConstruction.html_table_maker(sig_graphlets,"$(params["website"]["website_dir"])/_assets/$(params["website"]["page_name"])/sig_type_representations.html",imgs=sig_graphlets.Graphlet,figpath="../figs/")
        NetworkConstruction.html_table_maker(insig_graphlets,"$(params["website"]["website_dir"])/_assets/$(params["website"]["page_name"])/insig_type_representations.html",imgs=insig_graphlets.Graphlet,figpath="../figs/")  
    end 

    ##look at edge types in randomised networks
   # real_type_edgecounts = countmap(splat(tuple).(sort.(eachrow(hcat(map(x->vertexlist[x],first.(edgelist)),map(x->vertexlist[x],last.(edgelist)))))))
   # rand_types_edgecounts = map(y->(countmap(splat(tuple).(sort.(eachrow(hcat(map(x->y[x],first.(edgelist)),map(x->y[x],last.(edgelist)))))))),rand_types_set)
   # rand_edge_collection = vcat(collect.(rand_types_edgecounts)...)
   # rand_edge_df = DataFrame(graphlet = broadcast(first,rand_edge_collection),value = broadcast(last,rand_edge_collection))
   # random_edges = DataFrame()
   # for t in unique(rand_edge_df[!,:graphlet])
   #     rand_vals = filter(:graphlet=>x->x==t,rand_edge_df)[!,:value]
   #     rand_exp = sum(rand_vals)/N
   #     real_obs = real_type_edgecounts[t]
   #     append!(random_edges,DataFrame(Graphlet = first(t)*"_"*last(t)*"_edge", Expected = rand_exp,Observed = real_obs))       
   # end

    #pretty_table(random_edges,backend=:html,standalone = false)
    return hog_array
end 
export typed_representations

function load_relationships(file_path)  
    rel_dir = dirname(file_path)
    ## first divide csv into the node ids (integers) and the graphlet ids (string)
    #TODO adjustable to graphlet size?
    run(pipeline(`cut -d' ' -f 1-4 $(file_path)`,"$(rel_dir)/nodes"))
    run(pipeline(`cut -d' ' -f 5 $(file_path)`,"$(rel_dir)/graphlets"))
    
    ##First read in nodes
    node_path = rel_dir*"/nodes"
    graphlet_path = rel_dir*"/graphlets"
    @info "Loading relationship node info"
    nodes = CSV.File(node_path,delim=' ',header = false) |> Tables.matrix
    cleaner()
    @info "Loading relationship graphlet info"
    graphlets = CSV.File(graphlet_path,delim=' ',header = false) |> Tables.matrix
    cleaner()
    ##cleanup split files
    run(`rm $(node_path) $(graphlet_path)`)
    ##Getting by without splits for now, but maybe in future will be necessary
    ##see how many split files are required
    #bins = Int(ceil(filesize(node_path)/5368709120))
    ## find line splits to match bins (so we don't split across lines)
    #linecount = countlines(node_path) 
    #lines = Int(ceil(linecount/bins))
    ##if file is too big, we need to split up and load in separately
    #run(`split -l$(lines) $(node_path) $(rel_dir)/split_rels.`)  

    #nodes = [0 0 0 0]
    #for file in filter(x->occursin("split_rels.",x),readdir(dirname(file_path),join=true))
     #   Rel = CSV.read(file,DataFrame,delim = " ",header = false) 
    #cleaner()
     #  nodes = vcat(nodes,Matrix(Rel))
      ##cleanup split file
       # run(`rm $(file)`)
        #cleaner()
    #end
    #convert graphlets to vector instead of matrix
    graphlets = vec(graphlets)
    return [nodes, graphlets]
 end
export load_relationships

function kegg_information(network_counts)

    vertex_names = network_counts[!,:transcript_id]
    #get baseline entrez and kegg info about transcripts
    bio_dir = params["cache"]["bio_dir"]
    kegg_file = "$(bio_dir)/kegg_info.jld2"
    if (isfile(kegg_file))
        @info "Loading KEGG info from $kegg_file..."
        top_terms = cache_load(kegg_file,"top_terms")
    else
        @info "getting KEGG info..."
        top_terms = GraphletAnalysis.get_KEGG_pathways(vertex_names,"transcript")
        cache_save(kegg_file,["top_terms"=>top_terms ])
    end
    return top_terms

end
function go_information(network_counts)

    vertex_names = network_counts[!,:transcript_id]
    #get baseline GO info about transcripts
    bio_dir = params["cache"]["bio_dir"]
    go_file = "$(bio_dir)/go_info.jld2"
    if (isfile(go_file))
        @info "Loading GO info from $go_file..."

        top_terms = cache_load(go_file,"top_terms")
    else
        @info "getting GO info..."
        top_terms = GraphletAnalysis.get_GO_terms(vertex_names,"transcript")
        cache_save(go_file,["top_terms"=>top_terms ])
    end
    return top_terms

end

function ensembl_coverage(network_counts)
    names = network_counts.transcript_id
    bio_dir = params["cache"]["bio_dir"]
    entrez_file = "$(bio_dir)/entrez_map.jld2"
    if (isfile(entrez_file))

        entrez_map = cache_load(entrez_file,"entrez")
    else
        @info "getting entrez map..."
        entrez_map = GraphletAnalysis.get_entrez_ids(names,"transcript")
        cache_save(entrez_file,["entrez"=>entrez_map ])
    end
    return entrez_map
end

function biological_validation(network_counts)
    ##generate top terms for kegg and go
    ktt = kegg_information(network_counts)
    gtt = go_information(network_counts)
    ##at the moment return nothing. could return all components?
    return [ktt,gtt] 
end

function coincident_graphlets(network_counts,vertexlist,edgelist)
    #Coincident analysis
    coinc_dir = params["cache"]["coinc_dir"]
    run(`mkdir -p $(coinc_dir)`)
        

    #coincidents
    #TODO for smaller networks/lower order graphlets, it may take longer to load CSV than to just run
    #coincident analysis (set switch to turn on/off caching? detect automatically?)
    coincidents_file ="$coinc_dir/coincidents.csv"
    candidates_file = "$coinc_dir/candidates.jld2"
    if (isfile(coincidents_file) & isfile(candidates_file))
        @info "Loading coincidents dataframe from $coinc_dir..."
        Coincidents = CSV.read(coincidents_file,DataFrame)
        #because CSV converts the array columns to strings, we have to convert back (cost of using the easy/dirty CSV option!)
        fix(g) = split(replace(replace(replace(replace(g,("["=>"")),("]"=>"")),("\""=>"")),(" "=>"")),",")
        fix_int(g) = map(x->parse(Int,x),split(replace(replace(g,("["=>"")),("]"=>"")),","))
        fix_bool(g) = BitArray(map(x->parse(Int,x),split(replace(replace(replace(g,("["=>"")),("]"=>"")),("Bool"=>"")),",")))
        Coincidents.Vertices = fix_int.(Coincidents.Vertices)
        Coincidents.Entrez = fix_int.(Coincidents.Entrez)
        Coincidents.Ensembl = fix.(Coincidents.Ensembl)
        Coincidents.Ensembl = map(x->string.(x),Coincidents.Ensembl)
        Coincidents.Transcript_type = fix.(Coincidents.Transcript_type)
        Coincidents.Transcript_type = map(x->string.(x),Coincidents.Transcript_type)
        Coincidents.Inclusion = fix_bool.(Coincidents.Inclusion)

        ##load in matching candidate info
        @info "Loading candidate info from $coinc_dir..."
        entrez_id_vector = cache_load(candidates_file,"entrez_id_vector")
        candidates = cache_load(candidates_file,"candidates")
    else
        ##relationships 
        #Here we only load relationships file (which is huge!) if absolutely required i.e. the smaller resultant coincidents file does not exist
        graphlet_size = params["analysis"]["graphlet_size"]  
        rel_dir = params["cache"]["graphlet_enum_dir"]
        relationships_file = "$rel_dir/relationships.csv"
        if (isfile(relationships_file))
            @info "Loading per-edge graphlet relationships from $relationships_file"
            cleaner()
            rel,rel_types = load_relationships(relationships_file)
            cleaner()
        else
            #Setup workers
            distributed_setup(:ProjectFunctions,:GraphletCounting,:GraphletAnalysis,:NetworkConstruction)
            @info "Counting per-edge graphlet relationships..."
            ## graphlet_relationships already caches final file in csv form at temp_dir location
            Rel = GraphletCounting.graphlet_relationships(vertexlist,edgelist,graphlet_size,run_method="distributed",progress = true,temp_dir = rel_dir)
            ## remove workers now (no further distributed at this stage)
            rmprocs(workers())
            cleaner()
            rel,rel_types = load_relationships(relationships_file)
            cleaner()
            cleaner()
        end
        
        ##TODO remove need for 0 columns (each order is stored separately?) trim relationships list to just the graphlet order of interest. Zero column(s) will need to be trimmed for smaller graphlet orders atm
        
        

        rel = rel[findall(x->occursin(string(graphlet_size),x),rel_types),:]        
        rel_types = rel_types[findall(x->occursin(string(graphlet_size),x),rel_types)]
        if (size(rel)[2] > graphlet_size)
            rel = rel[:,(1+size(rel)[2]-graphlet_size):end]
        end
    
        ##get candidate information from kegg
        ##get baseline network info about biological function, with candidate info.
        @info "Getting candidate pathway info from KEGG..."
        vertex_names = network_counts[!,:transcript_id]
        entrez_id_vector, candidates,top_terms = GraphletAnalysis.get_KEGG_candidates(vertex_names,"transcript")   
        @info "Conducting per graphlet pathway coincidence analysis..."
        Coincidents = GraphletAnalysis.graphlet_coincidences(rel,rel_types,vertexlist,vertex_names,entrez_id_vector,candidates)
        @info "Saving coincidents and candidates at $coinc_dir..."

        CSV.write(coincidents_file,Coincidents)
        cache_save(candidates_file,["entrez_id_vector"=>entrez_id_vector,"candidates"=>candidates ])
    end
    
    ##get pathways as separate vector
    candidate_pathways = collect(keys(candidates))
    #orbit statistics
    #orbitsigs_file ="$val_dir/orbitsigs.jld"
    #if (isfile(orbitsigs_file))
    #        @info "Loading orbitsigs from $val_dir..."

    #collect dataframe for each node in this array
    #per_node_significance = Array{DataFrame,1}(undef,length(vertexlist))
    #choose just one order of graphlets (3 or 4)
    #TODO ideally getting sub_Coincidents shouldn't be necessary
    sub_Coincidents = filter(:Hom_graphlet=>x->occursin("$(params["analysis"]["graphlet_size"])-",x),Coincidents)


    orbit_dir = params["cache"]["orbit_dir"]
    run(`mkdir -p $(orbit_dir)`)
    orbit_sigs_file = "$orbit_dir/orbit_sigs.jld2" 
    ## select whether we are looking at "detailed" or "collated" significance
    orbit_sigs_method ="detailed"
    
    if (isfile(orbit_sigs_file))
        @info "Loading orbit significance dataframe from $orbit_dir..."
        orbit_sigs = cache_load(orbit_sigs_file,"orbit_sigs")
    else

        @info "Getting orbit significance statistics..."
        #table to showing whether each node (row) is included in each pathway (column)
        inkey = hcat([ in.(1:length(vertexlist),Ref(candidates[p])) for p in candidate_pathways ]...)
        if(orbit_sigs_method == "collated")
        end
        if(orbit_sigs_method == "detailed")
            orbit_sigs = @showprogress map(x->GraphletAnalysis.pernode_significance_detail(x,sub_Coincidents,params["analysis"]["graphlet_size"],candidate_pathways,inkey[x,:]),1:length(vertexlist))
        end
        cache_save(orbit_sigs_file,"orbit_sigs"=>orbit_sigs)
     end

    ## now compare the significance profile of those nodes that are not attached to a pathway to the average pathway profile of known pathway nodes
    ##convert to array form for comparisons
    orbit_sigs_array = map(x->Array(x[!,2:end]),orbit_sigs)
    #number of orbits
    last_col = size(orbit_sigs[1])[2]-1
    #define orbit names here for now
    orbit_names = names(orbit_sigs[1])[2:last_col+1]
    ##dictionary mapping number of times a node is in a pathway
    sig_pathway_occurences = countmap(vcat([candidates[x] for x in candidate_pathways]...))

    #ecdfs
    #store ecdf functions in table
    ecdf_table = Array{ECDF,2}(undef,length(candidate_pathways),last_col)
    for (i,p) in enumerate(candidate_pathways)
        for j in 1:last_col
            ecdf_table[i,j] = ecdf(map(x->x[i,j],orbit_sigs_array))
        end
    end
    #first check candidate probabilities across all categories
    known_pathway_probs = Array{Array{Float64,2},1}(undef,length(candidate_pathways))
    for (i,c) in enumerate(candidate_pathways)
        known_pathway_probs[i] = hcat([map(ecdf_table[i,j],map(x->x[i,j],orbit_sigs_array[candidates[c]])) for j in 1:last_col]...)
    end
    

    #collect known pathway vectors for corresponding known pathway nodes

    known_pathway_dfs = Array{DataFrame,1}(undef,length(candidate_pathways))
    for (i,p) in enumerate(candidate_pathways)
        subset = candidates[p]
        df_build = DataFrame(hcat(map(x->x[i,:],orbit_sigs_array[subset])...)',orbit_names)
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
    known_ecdf_table = Array{ECDF,2}(undef,length(candidate_pathways),last_col)
    for (i,p) in enumerate(candidate_pathways)
        known_ecdf_table[i,:] = ecdf.(eachcol(known_pathway_arrays[i]))
    end

    ##compare each unknown node to known ecdfs
    #set threshold for which to check on i.e. is node orbit count higher than the top ~thresh~% of known nodes
    thresh = 0.05
    ub = 1 -thresh
    #first check if there are any orbits for any pathways that have such low representation that 0 count exceed threshold. we will filter these out of analysis
    low_filter = (map(x->x(0),known_ecdf_table).<ub)
    #for each node, map each known ecdf to the corresponding orbit count and check against threshold, factoring in the low filter
    unknown_ecdf_comparison = map(y->(reshape(map(x->known_ecdf_table[x](y[x]),1:length(known_ecdf_table)),size(known_ecdf_table)[1],size(known_ecdf_table)[2]).>ub).*low_filter,orbit_sigs_array)
    #determine a node as significantly linked to a pathway if at least half its orbit counts are significant

    

    sig_cut = 1
    sig_check = map(x->(sum(x,dims=2).>(sig_cut)),unknown_ecdf_comparison)
    sig_nodes = findall(x->sum(x)>0,sig_check)
    sig_nodes_dict = Dict(Pair.(sig_nodes,map(x->candidate_pathways[vec(x)],sig_check[sig_nodes])))

    return sig_nodes_dict
end
export coincident_graphlets




