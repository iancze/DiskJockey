nchild = 4
keylist = Int[i for i=1:20]

# Chunk up the keylist into an acceptable number of keys per process.
max_keys_per_process = iceil(length(keylist)/nchild)
partial_proc = max_keys_per_process * nchild - length(keylist) # num of unfull processors
full_proc = nchild - partial_proc # number of full processors

# Feed max_keys_per_process keys to each of the full_proc, then feed
# (max_keys_per_process - 1) to each of the partial_proc
chunk_keylist = Array(Any, nchild)
start_key = 1
end_key = start_key + (max_keys_per_process - 1)
for i=1:full_proc
    println("Choosing from $start_key to $end_key")
    chunk_keylist[i] = keylist[start_key:end_key]
    start_key += max_keys_per_process
    end_key += max_keys_per_process
end

# Inch back to prep for filling the partially filled processes
end_key -= 1
# Now fill the partial_proc
for i=(full_proc+1):nchild
    println("Choosing from $start_key to $end_key")
    chunk_keylist[i] = keylist[start_key:end_key]
    start_key += (max_keys_per_process - 1)
    end_key += (max_keys_per_process - 1)
end

println(chunk_keylist)
