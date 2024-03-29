# SF-SFD: Stochastic training of Fourier coefficients for Space-Filling Designs

The method SF-SFD can be used to generate space-filling design at several common problem dimensions. It is used for tuning distribution functions in high-dimensional spaces in order to prevent concentration of measure for a finite sample. For more details about how to generate a sample, please read the section 'Basic Usage'. 

This technique directly addresses the issue of measure concentration and scales better to large dimensions than existing heuristic techniques such as Latin hypercube samples and low-discrepancy sequences such as the randomized Sobol sequence.

This work has been published in 2023 Winter Simulation Conference, if
you use it in your research, please consider citing our paper (see ``Citing``
section below)

## Dependencies

This package uses the latest versions of:
 - numpy (https://numpy.org) and
 - scipy (https://scipy.org)

Both can be installed using PyPI (``pip install`` command).

To run the tests, you will also need
 - pytest and
 - flake8

## Setup

Clone the main branch from GitHub

```
git clone https://github.com/sfdsampling/sfsfd
```

Then, either add the ``sfsfd`` directory to your Python packages, or
cd into the base directory and add it to your ``PYTHONPATH`` environment
variable:

```
cd sfsfd
export PYTHONPATH=$PYTHONPATH:`pwd`
```

## Testing

To check that it is installed, run our tests using

```
./run_tests.sh
```

or

```
flake8 sfsfd
pytest sfsfd
```

## Basic Usage

The key idea is to create an instance of the ``SamplingModel`` object,
then train a model using its ``initialize()`` method.

See the example script below:

```
""" Demonstrate the SamplingModel class on a 2D problem with 10 iterations
and a desired sample size of 30. """

from sfsfd import SamplingModel

dimension = 2
grid_cells_per_dimension = 3
iterations = 10
sample_size = 10
file_name = ('lhc2_dimension_' + str(dimension) + '_' + 'grid_cells_' +
             str(grid_cells_per_dimension)+'.txt')
model = SamplingModel( 
    dimension_of_input_space=dimension, 
    grid_size=grid_cells_per_dimension, 
    file_name=file_name,
    no_of_iterations_per_perturbation = iterations, # Increase with dimension
    sample_size = sample_size # This can be user's choice, e.g., 10
    )
# Train the model and save results to file_name
model.initialize()
```

## Experimental Results

To reproduce the results from our WSC 23 paper, see the ``README.md`` file
in the ``wsc_experiments`` subdirectory.

## Citing

To cite this work, please use the following

```
@inproceedings{garg2023sfsfd,
  author  = {Garg, Manisha and Chang, Tyler H. and Raghavan, Krishnan},
  title   = {{SF-SFD}: {S}tochastic optimization of {F}ourier coefficients for space-filling designs},
  booktitle = {Proc. 2023 Winter Simulation Conference (WSC 2023)},
  location = {San Antonio, TX, USA},
  note = {To appear},
  preprint = {https://arxiv.org/abs/2305.12043}
}
```

## Support

For questions, the primary author is:

 - Manisha Garg (manisha8@illinois.edu)
