push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/JudithExcalibur/")

@everywhere using parallel

## Master process

# Propose parameters, write down density, velocity, temperature structure.

## Workers

# Run RADMC3D, FFT, sample, evaluate lnprob, send back to master process

## Master process

# Decide parameters, advance chain

@everywhere function f(dset, p, key)
    0.0
    #println("Likelihood call to $dset, $p, $key")
end

pipes, dones = initialize(15, f)

println(procs())
println(workers())

distribute!(pipes, [0.0, 0.0])
result = gather!(pipes)

tic()
distribute!(pipes, [0.0, 0.0])

result = gather!(pipes)
toc()

println(result)

println("completing all results")

quit!(pipes)
#wait!(dones)

println("all results should be completed")

#wait!(dones)

println(procs())
#println(workers())



# # Initialization all children processes this way
# mypipe = Pipe(2)
# r2 = remotecall(2, brain, mypipe, 2)
#
# send!(mypipe, [1.0, 1.0, 1.2, 1.2])
# println("Sent message from master to child.")
#
# println("Master waiting for message retrieval from child")
# msg = get!(mypipe)
# println("Got $msg back!")
#
# send!(mypipe, "QUIT")
