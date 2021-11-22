
"wipe all buffers (this is a hard reset atm)
bufdo bwipeout!
" open/reload buffers to source and test .jl files
n src/*
n test/*
"this is homebase for now
b src/Construction


"set up local and remote terminals
term
file local
wincmd c
term zsh -is eval "~/git/rna-seq/nectar/launch-nectar-instance.sh graphlet-thread-test"
file remote
wincmd c

"setup unison file sync
term unison heterogeneous-graphlet-counting 
wincmd c

"open REPL terminal vertically on right
vert term zsh -is eval "conda activate nectar;~/git/rna-seq/nectar/launch-nectar-instance.sh graphlet-thread-test 'tmux a -t julia'"
file repl
wincmd p

"set up default window for vimslime to send to
let g:slime_default_config = {"bufnr": term_list()[0]}
let g:slime_dont_ask_default = 1
echo slime_default_config
echo term_list()
"for testing only
"e julia-remote.vim
