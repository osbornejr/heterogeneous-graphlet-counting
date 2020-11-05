mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
current_dir := $(notdir $(patsubst %/,%,$(dir $(mkfile_path))))
path_to_current_dir := $(patsubst %/,%,$(dir $(mkfile_path)))
julia_env := "osbornejr/julia-env:latest"

help: ## This help.
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

run: ## run julia
	docker run -it --rm -v $(path_to_current_dir):/home/osbornejr/app $(julia_env)  	

update: ## update conda and julia environments	
	docker image build --file docker/update/Dockerfile --tag $(julia_env) ./

build: ## build conda and julia environments from scratch
	docker image build --file docker/env/Dockerfile --tag $(julia_env) ./

push: ## push julia environment docker image
	docker push $(julia_env)

pull: ## pull julia environment docker image
	docker pull $(julia_env)
