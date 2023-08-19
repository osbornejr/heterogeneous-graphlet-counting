
function cache_remove(file::String)
    run(`rm -rf $(file)`)
end

function make_cache(dirs...;dir_name::String)
    #join args into dir path
    dir = join(dirs,"/") 
    #make dir if does not exist
    run(`mkdir -p $(dir)`)
    #store dir path in cache params
    params["cache"][dir_name] = dir
    ##remove contents if cache clear is set
    #TODO setup user input warning at start of run if any data is to be cleared.
    if(params["cache"]["clear"][dir_name])
        @info "Removing cached $(dir_name) data (and all dependent cached data) from $(dir)"
run(`rm -r $(dir)/`)
    end 
end

function cache_setup()
    ##create directory for run
    make_cache(dir_name="test_dir",cwd,params["cache"]["base_dir"],params["test_name"])
    ##switching to new method where each path is declared explicitly here (avoids run overlaps making a mess, and is easier with pluto load ins etc)
    ##data preprocessing dirs:
    make_cache(dir_name="round_dir",params["cache"]["test_dir"],"roundsig",string(params["data_preprocessing"]["roundsig"]))
    make_cache(dir_name="vst_dir",params["cache"]["round_dir"],"vst",string(params["data_preprocessing"]["vst"]))
    make_cache(dir_name="clean_dir",params["cache"]["vst_dir"],"expression_cutoff",string(params["data_preprocessing"]["expression_cutoff"]),"minreq",string(params["data_preprocessing"]["minreq"]),"clean_method",string(params["data_preprocessing"]["clean_method"]))
    make_cache(dir_name="norm_dir",params["cache"]["clean_dir"],"normalisation",params["data_preprocessing"]["norm_method"])
    make_cache(dir_name="sampling_dir",params["cache"]["norm_dir"],"sampling",string(params["data_preprocessing"]["variance_percent"]))
    #network construction dirs:
    ##at present we want to trial different bin sizes for the information measures, so we add the extra split here if "coexpression" is either MI or PID.
    if params["network_construction"]["coexpression"] in ["pidc","mutual_information"]
        
        make_cache(dir_name="similarity_dir",params["cache"]["sampling_dir"],"similarity",params["network_construction"]["coexpression"],params["network_construction"]["nbins"])
    else
        make_cache(dir_name="similarity_dir",params["cache"]["sampling_dir"],"similarity",params["network_construction"]["coexpression"])
    end
    

    make_cache(dir_name="adjacency_dir",params["cache"]["similarity_dir"],"threshold",string(params["network_construction"]["threshold"]),"threshold_method",params["network_construction"]["threshold_method"])

    #analyis dirs :
    if (params["network_construction"]["synthetic"] == true)
        make_cache(dir_name="anal_dir",params["cache"]["adjacency_dir"],"analysis","synthetic")
    else
        make_cache(dir_name="anal_dir",params["cache"]["adjacency_dir"],"analysis","real")
    end

    make_cache(dir_name="bio_dir",params["cache"]["anal_dir"],"biological_validation")
    make_cache(dir_name="community_dir",params["cache"]["anal_dir"],"communities")
    make_cache(dir_name="wgcna_dir",params["cache"]["community_dir"],"wgcna")
    make_cache(dir_name="graphlet_counting_dir",params["cache"]["anal_dir"],"graphlets",string(params["analysis"]["graphlet_size"]),"graphlet-counting")
    make_cache(dir_name="rep_dir",params["cache"]["graphlet_counting_dir"],"typed_representations","nullmodel",string(params["analysis"]["null_model_size"])*"_simulations")
    make_cache(dir_name="graphlet_enum_dir",params["cache"]["anal_dir"],"graphlets",string(params["analysis"]["graphlet_size"]),"graphlet-enumeration")
    make_cache(dir_name="coinc_dir",params["cache"]["graphlet_enum_dir"],"coincidents")
    make_cache(dir_name="orbit_dir",params["cache"]["coinc_dir"],"orbit-significance")
    #TODO setup archive method under new cache system
end

function cache_update(general::String,specific::Any)
    return cache_update(general,string(specific))
end
export cache_update

function cache_update(general::String,specific::String="",side_dir::String="")
    ##depreceated currently, as we prefer explicit directory creation at beginning of run
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
