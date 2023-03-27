#!/bin/bash
cd .. && export PYTHONPATH=$PYTHONPATH:`pwd` && cd wsc_experiments

# Run the tests with 10 random seeds ( 0--9 )
python3 compare_3d.py 0 1 2 3 4 5 6 7 8 9
python3 compare_4d.py 0 1 2 3 4
python3 compare_4d.py 5 6 7 8 9
python3 compare_5d.py 0 1
python3 compare_5d.py 2 3
python3 compare_5d.py 4 5
python3 compare_5d.py 6 7
python3 compare_5d.py 8 9
