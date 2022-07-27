#! /bin/bash

# The following README file can be run in the terminal to automatically create the selectivity plot
# In order to create the selectivity plot first set the paths for catmap then
# run a microkinetic model for all three cases of the H-CO barrier
for i in 0.0 0.1 -0.1
do
python3 scripts/run_catmap.py $i
done

# After the models ran figure 6b can be created running
python3 scripts/get_selectivity.py
