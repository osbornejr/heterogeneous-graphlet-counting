mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
current_dir := $(notdir $(patsubst %/,%,$(dir $(mkfile_path))))
path_to_current_dir := $(patsubst %/,%,$(dir $(mkfile_path)))
julia_env := "osbornejr/julia-env:latest"

help: ## This help.
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

run: ## run julia
	docker run --mount source=juliadotfolder,target=/home/osbornejr/.julia -it --rm -v $(path_to_current_dir):/home/osbornejr/app $(julia_env)  	

run-port: ## run julia with access to port 8000
	docker run -it --rm -p 8000:8000 -v $(path_to_current_dir):/home/osbornejr/app $(julia_env)  	

update: ## update conda and julia environments	
	docker image build --file docker/update/Dockerfile --tag $(julia_env) ./

build: ## build conda and julia environments from scratch
	docker image build --file docker/env/Dockerfile --tag $(julia_env) ./

full-build:## full build of all container dependencies. Needs to be done once in repository if the julia-env container is to be modified locally.
	docker image build --file docker/conda/Dockerfile --tag  conda-base:latest ./
	docker image build --file docker/julia/Dockerfile --tag julia:latest ./
	docker image build --file docker/env/Dockerfile --tag $(julia_env) ./

push: ## push julia environment docker image
	docker push $(julia_env)

pull: ## pull julia environment docker image
	docker pull $(julia_env)

http: ## create http server so that plots etc can be viewed on local machine
	python3 -m http.server 8787

sudo-docker: ##make docker usable as non-sudo (linux only)
	sudo groupadd docker
	sudo gpasswd -a $USER docker
	newgrp docker
