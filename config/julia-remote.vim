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

"set up local and remote terminals (without swapfiles)
nos term
file local
wincmd c
"nos term zsh -is eval "source ~/.conda/conda_init;conda activate nectar;~/git/rna-seq/nectar/launch-nectar-instance.sh graphlet-thread-test"
nos term ssh -t -i ~/.ssh/mit-derm joel@192.168.175.134 
file remote
wincmd c

"setup unison file sync
nos term ./bin/unison heterogeneous-graphlet-counting 
wincmd c

"setup terminal with port connection 
"nos term zsh -is eval "source ~/.conda/conda_init;conda activate nectar;~/git/rna-seq/nectar/launch-nectar-instance.sh graphlet-thread-test 'tmux a -d -t julia-server'"
nos term ssh -t -i ~/.ssh/mit-derm -L 8080:localhost:8000 joel@192.168.175.134 
file julia-port
wincmd c

"open REPL terminal vertically on right
""old nectar connecting command
"nos vert term zsh -is eval "source ~/.conda/conda_init;conda activate nectar;~/git/rna-seq/nectar/launch-nectar-instance.sh graphlet-thread-test 'tmux a -d -t julia'"
"nos vert term zsh -is eval "dtach -A /tmp/hgc-julia zsh -is eval 'make julia'"
nos vert term ssh -t -i ~/.ssh/mit-derm joel@192.168.175.134 
"commenting this out for now, cannot make it work smoothly
"'dtach -A /tmp/hgc-julia "cd /home/joel/git/heterogeneous-graphlet-counting && sudo make julia"'
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
