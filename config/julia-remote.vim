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

"set where remote is: nectar laptop race 
let remote_method="race"
let race_ip ="52-62-171-33"
let $RACE_ADDRESS = "ec2-user@ec2-" . race_ip . ".ap-southeast-2.compute.amazonaws.com"


"set up local and remote terminals (without swapfiles)
"""LOCAL terminal: to browse local files
nos term
file local
wincmd c

"""REMOTE terminal: to browse files on remote server
""nectar
"nos term zsh -is eval "source ~/.conda/conda_init;conda activate nectar;~/git/rna-seq/nectar/launch-nectar-instance.sh graphlet-thread-test"
""linux laptop
"nos term ssh -t -i ~/.ssh/mit-derm joel@192.168.175.134 
""race
execute 'nos term ssh -i ~/.ssh/race-hub ec2-user@ec2-' . race_ip . '.ap-southeast-2.compute.amazonaws.com'
file remote
wincmd c

"""UNISON terminal: to update files on remote server as they are edited locally
"setup unison file sync
nos term ./bin/unison heterogeneous-graphlet-counting-race
file unison
wincmd c

"""REMOTE PORT terminal: to connect with an open port (for pluto etc.)
"setup terminal with port connection 
""nectar
"nos term zsh -is eval "source ~/.conda/conda_init;conda activate nectar;~/git/rna-seq/nectar/launch-nectar-instance.sh graphlet-thread-test 'tmux a -d -t julia-server'"
""linux laptop
"nos term ssh -t -i ~/.ssh/mit-derm -L 8080:localhost:8000 joel@192.168.175.134 
""
""race
execute 'nos term ssh -i ~/.ssh/race-hub -L 8080:localhost:8000 ec2-user@ec2-' . race_ip . '.ap-southeast-2.compute.amazonaws.com'
file julia-port
wincmd c

"""REPL terminal: to run julia in on remote server
"open REPL terminal vertically on right
""nectar
"nos vert term zsh -is eval "source ~/.conda/conda_init;conda activate nectar;~/git/rna-seq/nectar/launch-nectar-instance.sh graphlet-thread-test 'tmux a -d -t julia'"
"nos vert term zsh -is eval "dtach -A /tmp/hgc-julia zsh -is eval 'make julia'"
""linux laptop
"nos vert term ssh -t -i ~/.ssh/mit-derm joel@192.168.175.134 
"commenting this out for now, cannot make it work smoothly
"'dtach -A /tmp/hgc-julia "cd /home/joel/git/heterogeneous-graphlet-counting && sudo make julia"'
""race
execute 'nos vert term ssh -i ~/.ssh/race-hub ec2-user@ec2-' . race_ip . '.ap-southeast-2.compute.amazonaws.com'
file repl
wincmd p

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
"set up default window for vimslime to send to (the last terminal created aka the repl
let g:slime_default_config = {"bufnr": term_list()[0]}
let g:slime_dont_ask_default = 1
"for testing only
"e julia-remote.vim


"add in infiltrator breakpoint here
nnoremap * O<Esc>xOMain.@infiltrate<Esc>:w<CR> 
"remove all breakpoints in file
nnoremap ** :%s/Main.@infiltrate//g<CR>:w<CR> 
