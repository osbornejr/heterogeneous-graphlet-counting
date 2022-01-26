import REPL
Base.atreplinit() do repl
	    # make the ; shell mode sticky
	         	!isdefined(repl, :interface) && (repl.interface = REPL.setup_interface(repl))
	             	repl.interface.modes[2].sticky=true
	          end
#using Weave
#using Debugger
using Revise
using Pkg
using Distributed
using ProjectFunctions
#Pkg.precompile()
#ENV["JULIA_PROJECT"] = "/home/osbornejr/app"
