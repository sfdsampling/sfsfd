from sfsfd import SamplingModel


def test_4d_sampler():
    """ Demonstrate the SamplingModel class on a 4D problem with 10 iterations
        and a desired sample size of 10. """

    dimension = 4
    grid_cells_per_dimension = 3
    iterations = 20
    sample_size = 10
    file_name = ('lhc2_dimension_' + str(dimension) + '_' + 'grid_cells_' +
                 str(grid_cells_per_dimension)+'.txt')
    model = SamplingModel( 
        dimension_of_input_space=dimension, 
        grid_size=grid_cells_per_dimension, 
        file_name=file_name,
        no_of_iterations_per_perturbation = iterations, # Increase with dim
        sample_size = sample_size # This can be user's choice, e.g., 10
        )
    # Train the model and save results to file_name
    model.initialize()


if __name__ == "__main__":
    " If called as main, run manually "

    test_4d_sampler()
