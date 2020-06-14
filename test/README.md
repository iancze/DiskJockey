# Testing Overview

The testing suites for DiskJockey are broken up into two separate areas. The tests contained within this main DiskJockey package are designed to be unit-like and cover the implemented modules and plotting scripts as efficiently as possible. They should be able to be run on a standard laptop in a meaningful amount of time. They are the continuous integration tests that will be run by Github Workflows.

You can run these tests by `cd`ing to the local repository of this repo and running

    julia> ]
    pkg> activate . 
    pkg> test

If you discover errors, please submit a [GitHub Issue](https://github.com/iancze/DiskJockey/issues).

## Output directory

There are a number of intermediate plots produced by the test suite, which are normally saved to a temporary directory on your filesystem that is automatically deleted when the code completes. If you want to save these files for later inspection, set the following environment variable in your shell 

    $ export JL_PLOT_TEMPDIR=your_directory

run the tests as normal, and a `your_directory` will be placed inside of the `tests` folder with the plot output from your tests.

# Cluster scale tests
That said, these tests are not sufficient to cover the functional integration required of running a full inference step. Since waiting for a failed launch on a cluster is a significant productivity hit, there is an additional suite of tests designed to test the actual inference loop. This does require a significant computational environment (since the parallelism is part of what needs to be tested), but the hope is that this can suite can run in a finite amount of time, and thus result in higher operational availability when run on actual datasets. These tests are available in the [DiskJockey-Tests.jl](https://github.com/iancze/DiskJockeyTests.jl) repository.