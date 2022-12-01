#!/bin/bash
export PYTHONPATH=$PYTHONPATH:`pwd`

# Run the tests
flake8 sfsfd
python3 sfsfd/tests/test_utils.py
python3 sfsfd/tests/test_2d_sample.py
python3 sfsfd/tests/test_4d_sample.py
