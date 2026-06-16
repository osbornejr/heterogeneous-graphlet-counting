# Heterogeneous Graphlet Counting on Coexpression Networks

This is a collection of Julia packages alongside a Docker and Conda environment which together can be deployed to identify and analyis graphlet counts on coexpression networks.

## Install
After cloning this repository (and all submodules), use the command

```
make full-build
```
to create the required Docker containers.

## Run

The code is best run from within a Julia REPL instance to allow for inspection of interim outputs.

Initiate the julia environment using the command
```
make julia
```


## Todo
- add stop/input for cache clearing option
- pluto cell organisation
- make script specific to publication results output
- bug in loading config file via params (depends on current working directory currently...)
- make config template for individual config files to check against
