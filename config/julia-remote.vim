"save current buffer number
let cur_buf = bufnr('%')

"wipe all buffers (this is a hard reset atm)
"bufdo bwipeout!
" open/reload buffers to source and test .jl files
"n src/*
"n test/*

"wipe any existing terminals before refreshing them
for i in term_list()
    exec "bd! ".i
endfor

"set where repl is: nectar laptop race local
let remote_method= "local"


if remote_method ==# "race"
    "for race, we need to set up depending on the ip of the currently used
    "instance.
    let race_ip ="13-236-165-237"
    let race_address = "ec2-user@ec2-" . race_ip . ".ap-southeast-2.compute.amazonaws.com"
    let $RACE_ADDRESS = race_address
    
    "unison in particular will need a separate config file depending on this
    "address. we generate that here first
    let remote_path = "/home/ec2-user/heterogeneous-graphlet-counting"
    let remote_host = race_address 
    let ssh_key = expand("~/.ssh/race-hub")
elseif remote_method ==# "laptop"
    "set unison config for linux laptop 
    let remote_path = "/home/joel/git/heterogeneous-graphlet-counting"
    let remote_host = "joel@192.168.175.134"
    let ssh_key = expand("~/.ssh/mit-derm")
elseif remote_method ==# "nectar"
    "set unison config for nectar 
    let remote_path = "/pvol/heterogeneous-graphlet-counting"
    let remote_host = "ubuntu@118.138.237.216"
    let ssh_key = expand("~/.ssh/nectar")
elseif remote_method ==# "local"
else
    echoerr "Unknown remote method: " . remote_method
    finish
endif

""do not need unison or remotes at all there is just a local repl
if remote_method ==# "local"
    "set up local and remote terminals (without swapfiles)
    """LOCAL terminal: to browse local files
    nos vert term
    file local
    wincmd p
else "set paths first for unison and remotes

    "set up local and remote terminals (without swapfiles)
    """LOCAL terminal: to browse local files
    nos term
    file local
    wincmd c

    let local_path = "/Users/osbornejr/git/heterogeneous-graphlet-counting"
    let remote_root = "ssh://" . remote_host . "//" . remote_path
    let config_path = local_path . "/config/unison/unison-" . remote_method . ".prf"
    
    "Build unison config based on method
    let base = readfile(expand("config/unison/heterogeneous-graphlet-counting.prf"))
    call insert(base, '#args for ssh', 0)
    call insert(base, 'sshargs = -i ' . ssh_key, 1)
    call insert(base, '#Roots', 2)
    call insert(base, 'root = ' . remote_root, 3)
    call insert(base, 'root = ' . local_path, 4)
    call insert(base, 'prefer = ' . local_path, 5)
    call insert(base, 'servercmd = ' . remote_path . '/bin/unison', 6)
    
    call writefile(base,config_path)
    
    "set unison config path to be where repo configs are stored
    let $UNISON = local_path . "/config/unison"
    
    
    """UNISON terminal: to update files on remote server as they are edited locally
    "setup unison file sync
    
    execute 'nos term ./bin/unison -ignorearchives unison-' . remote_method  
    file unison
    wincmd c
    
    
    "" for remote terminals this will depend on the remote_method
    if remote_method ==# "race"
        """REMOTE terminal: to browse files on remote server
        execute 'nos term ssh -i ~/.ssh/race-hub ec2-user@ec2-' . race_ip . '.ap-southeast-2.compute.amazonaws.com'
        file remote
        wincmd c
        
        
        """REMOTE PORT terminal: to connect with an open port (for pluto etc.)
        "setup terminal with port connection 
        execute 'nos term ssh -i ~/.ssh/race-hub -L 8080:localhost:8000 ec2-user@ec2-' . race_ip . '.ap-southeast-2.compute.amazonaws.com'
        file julia-port
        wincmd c
        
        """REPL terminal: to run julia in on remote server
        "open REPL terminal vertically on right
        execute 'nos vert term ssh -i ~/.ssh/race-hub ec2-user@ec2-' . race_ip . '.ap-southeast-2.compute.amazonaws.com'
        file repl
        wincmd p
    elseif remote_method ==# "laptop"
        """REMOTE terminal: to browse files on remote server
        nos term ssh -t -i ~/.ssh/mit-derm joel@192.168.175.134 
        file remote
        wincmd c
        
        
        """REMOTE PORT terminal: to connect with an open port (for pluto etc.)
        "setup terminal with port connection 
        nos term ssh -t -i ~/.ssh/mit-derm -L 8080:localhost:8000 joel@192.168.175.134 
        ""
        file julia-port
        wincmd c
        
        """REPL terminal: to run julia in on remote server
        "open REPL terminal vertically on right
        nos vert term ssh -t -i ~/.ssh/mit-derm joel@192.168.175.134 
        "commenting this out for now, cannot make it work smoothly
        "'dtach -A /tmp/hgc-julia "cd /home/joel/git/heterogeneous-graphlet-counting && sudo make julia"'
        file repl
        wincmd p
    elseif remote_method ==# "nectar"
        """REMOTE terminal: to browse files on remote server
        nos term zsh -is eval "source ~/.conda/conda_init;conda activate nectar;~/git/rna-seq/nectar/launch-nectar-instance.sh graphlet-thread-test"
        file remote
        wincmd c
        
        """REMOTE PORT terminal: to connect with an open port (for pluto etc.)
        "setup terminal with port connection 
        nos term zsh -is eval "source ~/.conda/conda_init;conda activate nectar;~/git/rna-seq/nectar/launch-nectar-instance.sh graphlet-thread-test 'tmux a -d -t julia-server'"
        file julia-port
        wincmd c
        
        """REPL terminal: to run julia in on remote server
        "open REPL terminal vertically on right
        nos vert term zsh -is eval "source ~/.conda/conda_init;conda activate nectar;~/git/rna-seq/nectar/launch-nectar-instance.sh graphlet-thread-test 'tmux a -d -t julia'"
        "nos vert term zsh -is eval "dtach -A /tmp/hgc-julia zsh -is eval 'make julia'"
        file repl
        wincmd p
    else
        echoerr "Unknown remote method: " . remote_method
        finish
    endif
    "end local or remote 
endif   

"clear old slime variable from each loaded buffer
let buffers = map(filter(copy(getbufinfo()), 'v:val.listed'), 'v:val.bufnr')
for i in buffers
    exec "b ".i
    if exists("b:slime_config")
        unlet b:slime_config
    endif
endfor

"this is homebase for now
exec ":b ".cur_buf
"set up default window for vimslime to send to (the last terminal created aka the repl if remote, or local if local.
let g:slime_default_config = {"bufnr": term_list()[0]}
let g:slime_dont_ask_default = 1
"for testing only
"e julia-remote.vim


"add in infiltrator breakpoint here
nnoremap * O<Esc>xOMain.@infiltrate<Esc>:w<CR> 
"remove all breakpoints in file
nnoremap ** :%s/Main.@infiltrate//g<CR>:w<CR> 
