import sampling_model
import numpy as np
def test_2d_sampler():
    """ Demonstrate the SamplingModel class on a 2D problem with 10 iterations
        and a desired sample size of 30. """

    dimension_list = np.arange(2,4,1)
    grid_cells_per_dimension = 3
    iterations = 10
    sample_size_list = np.arange(10,101,10)

    for dimension in dimension_list:
        for sample_size in sample_size_list:
            file_name = ('lhc_dimension_' + str(dimension) + '_' + 'grid_cells_' +
                        str(grid_cells_per_dimension)+'.txt')
            model = sampling_model.SamplingModel( 
                dimension_of_input_space=dimension, 
                grid_size=grid_cells_per_dimension, 
                file_name=file_name,
                no_of_iterations_per_perturbation = iterations, # Increase with dimension
                sample_size = sample_size # This can be user's choice, e.g., 10
                )
            # Train the model and save results to file_name
            model.initialize()


if __name__ == "__main__":
    " If run as main, call manually. "

    test_2d_sampler()
