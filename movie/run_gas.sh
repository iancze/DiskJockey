#!/bin/bash
#Launch all the relevant dust commands

for i in {0..9}
do julia make_gas.jl -r $i &
done
