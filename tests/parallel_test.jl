push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/JudithExcalibur/")
push!(LOAD_PATH, "/n/home07/iczekala/JudithExcalibur/")

println("On startup, we have ", nprocs(), " processes.")
println("On startup, we have ", nworkers(), " workers.")

@everywhere using parallel

@everywhere function f(dset, key, p)
    println("Likelihood call to $dset, $p, $key on process ", myid())
    return 0.0
end

@everywhere function initfunc(key::Int)
    println("Loading dset $key")
    return key
end

pipes = initialize(nchild, initfunc, f)

for j=1:nworkers()
    distribute!(pipes, [0.0, 0.0])
    result = gather!(pipes)
    println(result)
end

println("completing all results")
quit!(pipes)
println("all results should be completed")

# For some reason the workers() command is slow, and if you call it multiple times
# in quick succession it will give outdated information. Also, if you call it while
# there are no workers, it will return [1]
# for j=1:3
#     println(workers())
#     sleep(2)
# end
