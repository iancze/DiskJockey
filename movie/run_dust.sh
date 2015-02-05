#!/bin/bash
#Launch all the relevant dust commands

for i in {0..9}
do julia make_dust.jl -r $i &
done
