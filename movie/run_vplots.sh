#!/bin/bash
#Launch all the relevant dust commands

for i in {0..9}
do julia make_vis.jl -r $i > vis$i.out &
done
