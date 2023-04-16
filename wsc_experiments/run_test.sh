#!/bin/bash

cd .. && export PYTHONPATH=$PYTHONPATH:`pwd` && cd wsc_experiments

# Run small test probs with 1 random seed 0
python3 compare_dD.py 2   0
python3 compare_dD.py 3   0
python3 compare_dD.py 4   0
