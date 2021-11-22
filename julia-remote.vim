
"wipe all buffers (this is a hard reset atm)
"bufdo bwipeout!
" open/reload buffers to source and test .jl files
n src/*
n test/*


"wipe any existing terminals before refreshing them
for i in term_list()
    exec "bd! ".i
endfor

"set up local and remote terminals
term
file local
wincmd c
term zsh -is eval "conda activate nectar;~/git/rna-seq/nectar/launch-nectar-instance.sh graphlet-thread-test"
file remote
wincmd c

"setup unison file sync
term unison heterogeneous-graphlet-counting 
wincmd c

"open REPL terminal vertically on right
vert term zsh -is eval "conda activate nectar;~/git/rna-seq/nectar/launch-nectar-instance.sh graphlet-thread-test 'tmux a -t julia'"
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
b src/Construction
"set up default window for vimslime to send to (the last terminal created aka the repl
let g:slime_default_config = {"bufnr": term_list()[0]}
let g:slime_dont_ask_default = 1
"for testing only
"e julia-remote.vim
