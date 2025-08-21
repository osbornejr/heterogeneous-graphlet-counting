"""
    get_input_data()


Retrieve the raw input data associated with the curent experiment parameters.

## Example 
    julia> raw_counts = get_input_data()

"""
function get_input_data()
    file = join([params["cache"]["test_dir"],"raw_counts.jld2"],"/")
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

function get_input_data(config_file::String)
    load_config(config_file)
    get_input_data()
end

export get_input_data


"""
    get_output_data()


Retrieve all output data associated with the curent experiment parameters.

## Example 
    julia> raw_counts,processed_counts,similarity_matrix,adj_matrix,network_counts,vertexlist,edgelist = get_output_data()
"""
function get_output_data()
    ##method to get all outputs at once
    raw_counts =  get_input_data()
    samp_file = "$(params["cache"]["sampling_dir"])/sample_counts.jld2"
    sim_file = "$(params["cache"]["similarity_dir"])/similarity_matrix.jld2"
    if ((isfile(sim_file)) && (isfile(samp_file)))
        similarity_matrix = cache_load(sim_file,"similarity_matrix")
        processed_counts = cache_load(samp_file,"sample counts")
    else
        throw(ArgumentError("No cached files exist at either $sim_file or $samp_file, please run from scratch using run_all method."))
    end
    components,adj_matrix,network_counts,vertexlist,edgelist = get_network_construction()
    return [raw_counts,processed_counts,similarity_matrix,adj_matrix,network_counts,vertexlist,edgelist]
end

function get_output_data(config_file::String)
    load_config(config_file)
    get_output_data()
end

export get_output_data

"""
    get_preprocessed_data()


Retrieve all preprocessing data associated with the curent experiment parameters.

## Example
    julia> raw_counts,round_counts,vst_counts,clean_counts,norm_counts,sample_counts = get_preprocessed_data()
"""
function get_preprocessed_data()
    ##method to get preprocessed dataframes before any network construction
    raw_counts =  get_input_data()
    round_file = "$(params["cache"]["round_dir"])/round_counts.jld2"
    vst_file = "$(params["cache"]["vst_dir"])/vst_counts.jld2"
    clean_file = "$(params["cache"]["clean_dir"])/clean_counts.jld2"
    norm_file = "$(params["cache"]["norm_dir"])/norm_counts.jld2"
    samp_file = "$(params["cache"]["sampling_dir"])/sample_counts.jld2"
    if ( (isfile(round_file)) && (isfile(vst_file)) && (isfile(clean_file)) && (isfile(norm_file)) && (isfile(samp_file)) )
        round_counts = cache_load(round_file,"round counts")
        vst_counts = cache_load(vst_file,"vst counts")
        clean_counts = cache_load(clean_file,"clean counts")
        norm_counts = cache_load(norm_file,"norm counts")
        sample_counts = cache_load(samp_file,"sample counts")
    else
        throw(ArgumentError("No cached files exist at least one of $round_file, $vst_file, $clean_file, $norm_file or $samp_file, please run from scratch using run_all method."))
    end
    return [raw_counts,round_counts,vst_counts,clean_counts,norm_counts,sample_counts]
end

function get_preprocessed_data(config_file::String)
    load_config(config_file)
    get_preprocessed_data()
end

export get_preprocessed_data


"""
    get_similarity_matrix()


Retrieve the similarity matrix associated with the curent experiment parameters.

## Example 
    julia> similarity_matrix = get_similarity_matrix()
"""
function get_similarity_matrix() 
    ##alt method that allows loading of cache if output exists, and gives an error otherwise.
    sim_file = "$(params["cache"]["similarity_dir"])/similarity_matrix.jld2"
    if (isfile(sim_file)) 
        sim_matrix = cache_load(sim_file,"similarity_matrix")
   else
       throw(ArgumentError("No cached file exists at $sim_file, please run experiment to generate similarity matrix."))
    end

    return sim_matrix      
end

function get_similarity_matrix(config_file::String)
    load_config(config_file)
    get_similarity_matrix()
end

export get_similarity_matrix

"""
    get_network_construction()


Retrieve all network information associated with the curent experiment parameters.

## Example 
    julia> components,adj_matrix,network_counts,vertexlist,edgelist = get_network_construction()
"""
function get_network_construction()
    ##alt method that allows loading of cache if output exists, and gives an error otherwise.
    samp_file = "$(params["cache"]["sampling_dir"])/sample_counts.jld2"
    adj_file = "$(params["cache"]["adjacency_dir"])/adjacency_matrix.jld2"
    if ((isfile(adj_file)) && (isfile(samp_file)))
        pre_adj_matrix = cache_load(adj_file,"pre-adj_matrix")
        adj_matrix = cache_load(adj_file,"adjacency_matrix")
        components = cache_load(adj_file,"components")
        sample_counts = cache_load(samp_file,"sample counts")
        # need to get includes again here based on select_components option
        includes = select_network_components(components)
   else
       throw(ArgumentError("No cached files exist at either $adj_file or $samp_file, please provide processed counts input data."))
    end

    network_counts = sample_counts[includes,:]    
    #maintain list of vertices in graph
    vertexlist = copy(network_counts[!,:transcript_type])     
    edgelist = NetworkConstruction.edgelist_from_adj(adj_matrix)

    #recalculate component labels here so that they match network indices, not processed count indices (at present this is the only place where components are exposed, so can happen here even if cached components are indexed on processed counts)
    updated_components = NetworkConstruction.network_components(adj_matrix)
    return updated_components,adj_matrix,network_counts,vertexlist,edgelist       
end

function get_network_construction(config_file::String)
    load_config(config_file)
    get_network_construction()
end

export get_network_construction

"""
    get_communities()
Retrieve community information associated with the curent experiment parameters.

## Example 
    julia> communities = get_communities()
"""
function get_communities()
    communities_file ="$(params["cache"]["community_dir"])/communities.jld2"

    if (isfile(communities_file))
        comms = cache_load(communities_file,"communities")
    else
       throw(ArgumentError("No cached file exists at $communities_file, please run community analysis."))
    end
    return comms
end

function get_communities(config_file::String)
    load_config(config_file)
    get_communities()
end

export get_communities

"""
    get_wgcna()

Retrieve WGCNA network and community information associated with the curent experiment parameters.

## Example 
    julia> wgcna_network,wgcna_comm = get_wgcna()
"""
function get_wgcna()

    wgcna_file = "$(params["cache"]["wgcna_dir"])/wgcna.jld2"
    if (isfile(wgcna_file))
        wgcna_network = cache_load(wgcna_file,"adjacency")
        wgcna_comm = cache_load(wgcna_file,"communities")
    else
       throw(ArgumentError("No cached file exists at $wgcna_file, please run wgcna analysis."))
    end
    return [wgcna_network,wgcna_comm]
end

function get_wgcna(config_file::String)
    load_config(config_file)
    get_wgcna()
end

export get_wgcna

"""
    get_biological_validation()

Retrieve basic biological validation information associated with the curent experiment parameters.

## Example 
    julia> kegg_top_terms,go_top_terms = get_biological_validation()
 """
function get_biological_validation()
    #get baseline entrez and kegg info about transcripts
    bio_dir = params["cache"]["bio_dir"]
    kegg_file = "$(bio_dir)/kegg_info.jld2"
    if (isfile(kegg_file))
        ktt = cache_load(kegg_file,"top_terms")
    else
            throw(ArgumentError("No cached file exists at $kegg_file, please ensure biological validation has been run on network."))
    end
    go_file = "$(bio_dir)/go_info.jld2"
    if (isfile(go_file))
        gtt = cache_load(go_file,"top_terms")
    else
            throw(ArgumentError("No cached file exists at $go_file, please ensure biological validation has been run on network."))
    end
    return [ktt,gtt]
end

function get_biological_validation(config_file::String)
    load_config(config_file)
    get_biological_validation()
end

export get_biological_validation

"""
    get_graphlet_counts()

Retrieve graphlet counts associated with the network under curent experiment parameters.
Also includes time that run took for downstream analysis.

## Example
    julia> graphlet_counts,timer = get_graphlet_counts()
 """
function get_graphlet_counts()

    anal_dir = params["cache"]["graphlet_counting_dir"]    
    graphlet_file = "$anal_dir/graphlet-counting.jld2" 
    if (isfile(graphlet_file))
        graphlet_counts = cache_load(graphlet_file,"graphlets")
        timer = cache_load(graphlet_file,"time")
    else
            throw(ArgumentError("No cached file exists at $graphlet_file, please ensure graphlet counting has been run on network."))
    end
    return [graphlet_counts,timer]
end

function get_graphlet_counts(config_file::String)
    load_config(config_file)
    get_graphlet_counts()
end

export get_graphlet_counts
