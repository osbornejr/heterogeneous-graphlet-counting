mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
current_dir := $(notdir $(patsubst %/,%,$(dir $(mkfile_path))))
path_to_current_dir := $(patsubst %/,%,$(dir $(mkfile_path)))

help: ## This help.
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

run: ## run julia
	docker run -it --rm -v $(path_to_current_dir):/home/osbornejr/app julia-env:latest  	

update: ## update conda and julia environments	
	docker image build --file docker/update/Dockerfile --tag julia-env:latest ./

build: ## build conda and julia environments from scratch
	docker image build --file docker/env/Dockerfile --tag julia-env:latest ./
