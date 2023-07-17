# Experiments

This subdirectory contains information and scripts for replicating
results in our paper

*Garg et al.
SF-SFD: Stochastic optimization of Fourier coefficients to generate
space-filling designs.
In WSC 2023.*

## Critical Scripts

To reproduce our experiments, the ``compare_dD.py`` script can be
called with a problem dimension followed by a list of random seeds.
It will perform the experiment for all sample sizes in the range
``n = 100, 200, ..., 500``.

Before running, make sure that SF-SFD has been added to the ``PYTHONPATH``,
following the instructions in the top-level ``README.md``.

For example, in the bash shell, to perform experiments at dimension 2 with
random seeds 11, 22, and 33, use the following command:

```
python compare_dD.py 2 11 22 33
```

To quickly test that you have dependencies installed, use
```
./run_test.sh
```

To reproduce our experiments, use
```
./run_all.sh
```

## Subdirectory Structure

 - The ``results`` directory contains the raw data that we collected
   when running our experiments for the paper.
 - The ``plots`` directory contains two scripts for reproducing the
   plots and tables in Section 5, based on the raw data in ``results``.

