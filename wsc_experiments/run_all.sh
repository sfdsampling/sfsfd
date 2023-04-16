#!/bin/bash

cd .. && export PYTHONPATH=$PYTHONPATH:`pwd` && cd wsc_experiments

# Run the scaling tests with 10 random seeds ( 0--9 )
python3 compare_dD.py 5   0 1 2 3 4 5 6 7 8 9
python3 compare_dD.py 10  0 1 2 3 4 5 6 7 8 9
python3 compare_dD.py 15  0 1 2 3 4 5 6 7 8 9
python3 compare_dD.py 20  0 1 2 3 4 5 6 7 8 9
python3 compare_dD.py 25  0 1 2 3 4 5 6 7 8 9
python3 compare_dD.py 30  0 1 2 3 4 5 6 7 8 9
