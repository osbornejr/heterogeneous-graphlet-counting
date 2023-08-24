"""
    get_input_data()


Retrieve the raw input data associated with the curent experiment parameters.

Output is in form: `raw_counts`

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
export get_input_data


"""
    get_output_data()


Retrieve all output data associated with the curent experiment parameters.

Output is in form: 

`raw_counts`,`processed_counts`,`similarity_matrix`,`adj_matrix`,`network_counts`,`vertexlist`,`edgelist`
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
export get_output_data

"""
    get_preprocessed_data()


Retrieve all preprocessing data associated with the curent experiment parameters.

Output is in form: `raw_counts`,`round_counts`,`vst_counts`,`clean_counts`,`norm_counts`,`sample_counts`
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
export get_preprocessed_data

"""
    get_network_construction()


Retrieve all network information associated with the curent experiment parameters.

Output is in form: `components`,`adj_matrix`,`network_counts`,`vertexlist`,`edgelist`
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
   else
       throw(ArgumentError("No cached files exist at either $adj_file or $samp_file, please provide processed counts input data."))
    end
    #Trim nodes with degree zero

    largest = findmax(length.(components))[2]
    network_counts = sample_counts[components[largest],:]    
    #maintain list of vertices in graph
    vertexlist = copy(network_counts[!,:transcript_type])     
    edgelist = NetworkConstruction.edgelist_from_adj(adj_matrix)
    return components,adj_matrix,network_counts,vertexlist,edgelist       
end
export get_network_construction


"""
    get_wgcna()

Retrieve WGCNA network and community information associated with the curent experiment parameters.

Output is in form: `wgcna_network`,`wgcna_comm`
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
export get_wgcna

"""
    get_biological_validation()

Retrieve basic biological validation information associated with the curent experiment parameters.

Output is in form: `kegg_top_terms`,`go_top_terms`
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
export get_biological_validation
