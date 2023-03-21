from sfsfd import sampling_model
import logging
import numpy as np
from scipy.stats import qmc
from scipy.spatial.distance import pdist

# Problem hyperparams
dimension_list = np.arange(4,5,1) # 1 levels
grid_cells_per_dimension = 3 # discretization level of 3
sample_size_list = np.arange(10,101,30) # 4 levels
no_of_iterations_per_perturbation = 10 # starting number of perturbations
adaptive_sample_size = 50 # increase by 1 every 100 iterations

def comparison():
    """ Main driver routine that performs comparison of all 3 techniques. """

    # Activate info-level logging
    logging.basicConfig(level=logging.INFO,
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')

    np.random.seed(31423) # Set numpy random seed (based on date 3.14.23)
    # Loop over all valid dimensions and sample sizes and generate data
    for dimension in dimension_list:
        for sample_size in sample_size_list:
            file_name = ('comparison_dimension_' + str(dimension)+'.txt')
            with open(file_name, "a") as file_instance:
                #file_instance.write("The methods compared are: Latin Hypercube," +
                #            " Sobol sequences and "+
                #            " Sf-sfd\n"
                #            )
    
                file_instance.write("********************************\n")
                file_instance.write("           SF-SFD                \n")
                file_instance.write("********************************\n")

            sfd_sample(dimension, sample_size, file_name)
            
            latin_hypercube(dimension, sample_size, file_name)
            
            sobol_seq(dimension, sample_size, file_name)


def sfd_sample(dimension, sample_size, file_name):
    """ Demonstrate the SamplingModel class on a 2D problem with 10 iterations
        and a desired sample size of 30. """

    grid_cells_per_dimension = 3
    file_name = file_name
    model = sampling_model.SamplingModel( 
        dimension_of_input_space=dimension, 
        grid_size=grid_cells_per_dimension, 
        file_name=file_name,
        no_of_iterations_per_perturbation = no_of_iterations_per_perturbation, # Increase with dimension
        adaptive_sample_size = adaptive_sample_size,
        sample_size = sample_size # This can be user's choice, e.g., 10
        )
    # Train the model and save results to file_name
    model.initialize()

### Other sampling methods from scipy.qmc ###
              
def latin_hypercube(dimension, sample_size, file_name):
    sampler = qmc.LatinHypercube(dimension)
    sample = sampler.random(sample_size)
            
    dis_lhs = qmc.discrepancy(sample)
    maximin_lhs = maximindist(sample)
    e_optimality_lhs = e_optimality(sample)
    weighted_criteria_lhs = dis_lhs - maximin_lhs - e_optimality_lhs
    with open(file_name, "a") as file_instance:
        file_instance.write("********************************\n")
        file_instance.write("           LHS                  \n")
        file_instance.write("********************************\n")

        #file_instance.write("The sample is = " +
        #                        f"{sample}\n")
        file_instance.write("The discrepancy of the sample is =" +
                                f"{dis_lhs}\n")
        file_instance.write("The maximin distance of the sample is =" +
                                f"{maximin_lhs}\n")
        file_instance.write("The min eigenvalue of the sample is =" +
                                f"{-e_optimality_lhs}\n")
        file_instance.write("The weighted criterion value is =" +
                                f"{weighted_criteria_lhs/3}\n")

def sobol_seq(dimension, sample_size, file_name):
    sampler = qmc.Sobol(dimension, scramble=False)
    sample = sampler.random(sample_size)
    dis_sobol = qmc.discrepancy(sample)
    maximin_sobol = maximindist(sample)
    e_optimality_sobol = e_optimality(sample)
    weighted_criteria_sobol = dis_sobol-maximin_sobol-e_optimality_sobol
    with open(file_name, "a") as file_instance:
        file_instance.write("********************************\n")
        file_instance.write("           SOBOL                 \n")
        file_instance.write("********************************\n")

        #file_instance.write("The sample is = " +
        #                        f"{sample}\n")
        file_instance.write("The discrepancy of the sample is =" +
                                f"{dis_sobol}\n")
        file_instance.write("The maximin distance of the sample is =" +
                                f"{maximin_sobol}\n")
        file_instance.write("The min eigenvalue of the sample is =" +
                                f"{-e_optimality_sobol}\n")
        file_instance.write("The weighted criterion value is =" +
                                f"{weighted_criteria_sobol/3}\n")

### Helper function to calculate criteria scores ###

def maximindist(sample):
    maximindistance = min(pdist(sample)) # By default Euclidean distance
    return maximindistance

def e_optimality(sample):
    sample_arr = np.array([arr.tolist() for arr in sample])
    t = sample_arr.T # 4x10
    u,s,v = np.linalg.svd(t)
    min_eigenvalue = np.min(s)
    return min_eigenvalue

### Driver call ###

if __name__ == "__main__":
    " If run as main, call driver. "

    comparison()

