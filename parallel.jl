module parallel

export Pipe, send!, get!, brain, initialize, distribute!, gather!, quit!

# Use this code in a separate module by doing
# @everywhere using parallel

# Pipe is designed to act like Python's two-ended PIPE
# which means master and child can both put and take things off the queue
type Pipe
    m::Int #Process ID of master
    c::Int #Process ID of child
    mtc::RemoteRef # Used to send message from Master to Child
    ctm::RemoteRef # Used to send message from Child to Master
end

# Automatically create the remote refs by simply specifying the child process ID
function Pipe(child::Int, master::Int = 1)
    return Pipe(master, child, RemoteRef(master), RemoteRef(child))
end

# Send a message to the pipe. Automatically figures out which RemoteRef to use
# based on whether the master or child calls it
function send!(pipe::Pipe, msg)
    id = myid()
    if id == pipe.m
        put!(pipe.mtc, msg)
    elseif id == pipe.c
        put!(pipe.ctm, msg)
    else
        println("Wrong process communication. ID's current: ", id, " master: ", pipe.m,
        " child: ", pipe.c)
    end
end

# Take a message from the pipe. Automatically figures out which RemoteRef to use
# based on whether the master or child calls it
function get!(pipe::Pipe)
    id = myid()
    if id == pipe.m
        return take!(pipe.ctm)
    elseif id == pipe.c
        return take!(pipe.mtc)
    else
        println("Wrong process communication. ID's current: ", id, " master: ", pipe.m,
        " child: ", pipe.c)
    end
end

# This function is started on each child process, looping constantly. It waits
# until the pipe has put any new information on (ie, the proposed parameters).
# Once information appears, the child takes the value, does its work and then
# puts the answer back on the pipe for the main process to consume. Then it goes
# back to waiting for more information.
# initfunc is designed to return the dataset, or any object that gets passed
# to f, specific to this process
function brain(pipe::Pipe, key::Int, initfunc::Function, f::Function)
    id = myid()
    println("Initialized brain $id with $key")

    #Load the dataset according to this key
    dset = initfunc(key)

    while true
        p = get!(pipe)
        #println("Process $id received parameters $p")

        #Check to see if the message is actually the quit signal
        if p == "QUIT"
            println("Process $id completing!")
            return

        # otherwise proceed with the likelihood calculation
        else
            # Likelihood function is declared elsewhere and takes as input the dataset,
            # the parameters, and the key
            ans = f(dset, key, p)

            # Return the answer to the main process
            send!(pipe, ans)
        end
    end
end

# All of the following functions are run ONLY the master process

# Set up nchild child processes and pipes and return an array of the pipes
# function takes three arguments: the dataset, integer key, and vector of parameters
function initialize(nchild::Int, initfunc::Function, f::Function)

    pipes = Array(Pipe, nchild)

    # The master process is always labelled 1, so nprocessors = nchild + 1
    # yet the pipes array is indexed from 1 to nchild, hence the offset
    for i=1:nchild
        child_id = i + 1
        pipes[i] = Pipe(child_id) # Create a pipe from the master process to this child ID

        # Initialize the process with the key, init function, and the function
        remotecall(child_id, brain, pipes[i], child_id, initfunc, f)
    end

    return pipes
end

# Distribute parameters to all sub processes for likelihood evaluation
function distribute!(pipes::Vector{Pipe}, p::Vector{Float64})
    for pipe in pipes
        send!(pipe, p)
    end
end

# Collect the results from all of the sub processes and sum together
function gather!(pipes::Vector{Pipe})
    return sum([get!(pipe) for pipe in pipes])
end

# Call this at the very end of the run to close all the pipes
function quit!(pipes::Vector{Pipe})
    for pipe in pipes
        send!(pipe, "QUIT")
    end
    #Kill all the processes
    rmprocs(workers())
end

# might be worthwhile to investigate the driver code mentioned in the manual. Seems like more of an MPI startup approach (everything is loaded everywhere, and we can check to see if we are master or not).


end #module
