#!/bin/bash
cd .. && export PYTHONPATH=$PYTHONPATH:`pwd` && cd wsc_experiments

# Run the tests
python3 compare_3d.py
python3 compare_4d.py
python3 compare_5d.py
