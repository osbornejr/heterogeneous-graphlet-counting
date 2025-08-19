using YAML

function load_config(config_file::String;template::String = "config/run-files/template.yaml",no_template = false)
    if no_template
        ##just take input from supplied file directly (this is the old method)
        #setup parameters to be read.writable for all functions in this module
        global params = YAML.load_file(config_file)
        cache_setup()
    else
        #first use template
        param_file = YAML.load_file(template)
        #load in supplied config
        supplied_file = YAML.load_file(config_file)
        #from supplied, get test_name that corresponds to experiment of interest
        experiment = supplied_file["test_name"]
        #load experiment specific parameters in case some aren't in supplied
        experiment_file = YAML.load_file("config/run-files/$(experiment)/$(experiment).yaml")
        #merge parameters from supplied and experiment with config
        #first we do experiment file
        recursive_merge!(param_file,experiment_file)
        #then override with any run specific options
        recursive_merge!(param_file,supplied_file)

        
        #setup parameters to be read.writable for all functions in this module
        global params = param_file
        cache_setup()
    end
end
export load_config
export params

function recursive_merge!(base::Dict, overrides::Dict)
    for (k, v) in overrides
        if haskey(base, k) && isa(v, Dict) && isa(base[k], Dict)
            recursive_merge!(base[k], v)
        else
            base[k] = v
        end
    end
    return base
end
