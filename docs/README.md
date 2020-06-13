# Building the docs

As of Julia v1.4, this is the following plan to build the docs.

    # install the necessary packages to the docs environment
    julia> ]
    pkg> activate .
    pkg> build 

    # run the make.jl script assuming this environment 
    # --project assumes Project.toml of current directory
    $ julia --project make.jl

Output will be in `docs/build/`