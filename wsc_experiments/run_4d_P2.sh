#!/bin/bash

#COBALT -n 1
#COBALT -q gpu_smx2_v100
#COBALT -A Performance
#COBALT -t 6:00:00
#COBALT -I

source activate sfd

export http_proxy="http://0.0.0.0:8889"
export https_proxy="https://0.0.0.0:8889"
export ftp_proxy="ftp://0.0.0.0:8889"



cd .. && export PYTHONPATH=$PYTHONPATH:`pwd` && cd wsc_experiments

# Run the tests with 10 random seeds ( 0--9 )
python3 compare_4d.py 5 6 7 8 9