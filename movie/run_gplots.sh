#!/bin/bash
#Launch all the relevant dust commands

for i in {0..9}
do julia make_gplots.jl -r $i > gas$i.out &
done
