mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
current_dir := $(notdir $(patsubst %/,%,$(dir $(mkfile_path))))
path_to_current_dir := $(patsubst %/,%,$(dir $(mkfile_path)))
julia_env := "osbornejr/julia-env:latest"
uname := $(shell uname)

help: ## This help.
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

julia: ## run julia
	docker run --mount source=dot-julia,target=/home/osbornejr/.julia -e startup=config/startup-files/env-startup.jl -it --rm -u $(id -u):$(id -g) -v $(path_to_current_dir):/home/osbornejr/app $(julia_env)  	

julia-clean: ## run julia with no startup 
	docker run --mount source=dot-julia,target=/home/osbornejr/.julia -it --rm  -u $(id -u):$(id -g) -v $(path_to_current_dir):/home/osbornejr/app $(julia_env)  	

julia-port: ## run julia with access to port 8000
	docker run --mount source=dot-julia,target=/home/osbornejr/.julia -it --rm -u $(id -u):$(id -g) -p 8000:8000 -v $(path_to_current_dir):/home/osbornejr/app $(julia_env)  	

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
	#sudo groupadd docker
	sudo gpasswd -a $(USER) docker
	newgrp docker
UNAME:= $(shell uname)
ifeq ($(UNAME),Darwin)
unison_file := "unison-v2.51.4+ocaml-4.05.0+x86_64.macos-10.15.tar.gz"
endif
ifeq ($(UNAME),Linux)
unison_file := "unison-v2.51.4+ocaml-4.05.0+x86_64.linux.tar.gz"
endif
unison: ##use this to sync repo with a remote host. (this command just installs unison)
	rm -rf bin/unison
	mkdir -p bin/unison
	##download file
	wget --no-check-certificate --content-disposition -P ./bin/unison/ "https://github.com/bcpierce00/unison/releases/download/v2.51.4/$(unison_file)"
	##extract and tidy
	tar -xvzf bin/unison/$(unison_file) -C bin/unison
	mkdir bin/temp
	mv bin/unison/bin/unison* bin/temp/
	rm -r bin/unison
	mv bin/temp/unison* bin/
	rm -r bin/temp
	#make sure config file is in right place
	#mkdir -p $HOME/.unison
	#cp config/unison/heterogeneous-graphlet-counting.prf $HOME/.unison/
	#set up share directory
	mkdir -p output/share
	# if on mac, set up fsmonitor
	@if [ "$$(uname)" = "Darwin" ]; then \
		brew install autozimu/homebrew-formulas/unison-fsmonitor; \
	fi
nectar-connect: ## ssh into nectar VM   


swapfile: ##manually create large (5g) swapfile for kernel to use
	###DO NOT RUN IF SWAPFILE ALREADY EXISTS!! TODO add check for this
	## create swapfile
	sudo dd if=/dev/zero of=/mnt/swapfile bs=1M count=5120
	#format
	sudo mkswap /mnt/swapfile
	#add to swap
	sudo swapon /mnt/swapfile

edit: ## vim command to initialise editing environment via nectar-ssh
	vim -c "source config/julia-remote.vim" packages/ProjectFunctions/src/* packages/GraphletCounting/src/* packages/GraphletAnalysis/src/* packages/NetworkConstruction/src/* packages/DataPreprocessing/src/* config/run-files/*

website: ##set up website branch under separate subdir
	mkdir -p website
	git fetch origin website
	git worktree website website


tmux_session:="julia-vim-ssh"
tmux: ## setup tmux session for this project based on edit command
	tmux new -d -s $(tmux_session)
	tmux send-keys -t $(tmux_session) 'make edit' C-m C-m
	tmux set-option -t $(tmux_session) status-bg blue
	tmux set-option -t $(tmux_session) status-left '#S'
	tmux set-option -t $(tmux_session) status-left-length 100
	tmux set-window-option -t $(tmux_session) window-status-current-format ''
