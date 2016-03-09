# Priors

When fitting certain disks, it may be worthwhile to include new information from separate analyses. For example, when fitting multiple CO transitions, it may be worthwhile to make priors that constrain the temperature profile.

Because flexible priors like these are difficult to "hard-code" into the package, there is an additional functionality that allows the user to write their own prior, which will overwrite the default prior in the package.

You can copy a sample stub to your current working directory via

    $ DJInitialize.jl --prior

Then, open this file with your favorite text editor. It is important that you mimic the exact same function call as in `src/model.jl`.
