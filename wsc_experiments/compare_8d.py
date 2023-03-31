from sfsfd import sampling_model
import logging
import numpy as np
from scipy.stats import qmc
from scipy.spatial.distance import pdist
import csv

# Problem hyperparams
dimension_list = np.arange(8,9,1) # 1 levels
grid_cells_per_dimension = 10 # discretization level of 10
sample_size_list = np.arange(10,101,30) # 4 levels
no_of_iterations_per_perturbation = 25 # starting number of perturbations
adaptive_sample_size = 25 # increase by 1 every 25 iterations
fieldnames = ["method","discrepancy","maximin",
              "eigen_value","cumulative", "time_to_sol"]
weights = np.ones(3)
weights[1] = 0.1
weights[2] = 0.01
weights = weights / np.sum(weights)

def comparison(iseed, file_name_csv):
    """ Main driver routine that performs comparison of all 3 techniques. """

    # Activate info-level logging
    logging.basicConfig(level=logging.INFO,
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')

    # Loop over all valid dimensions and sample sizes and generate data
    for dimension in dimension_list:
        with open(file_name_csv, "a") as file_instance:
            file_instance.write(f"Dimension: {dimension}\n")
            file_instance.write(f"\nSample seed: {iseed}\n\n")
            writer = csv.DictWriter(file_instance, fieldnames=fieldnames)
            writer.writeheader()
        for sample_size in sample_size_list:
            with open(file_name_csv, "a") as file_instance:
                file_instance.write(f"\nSample size: {sample_size}\n\n")
        
            np.random.seed(iseed) # Set numpy random seed
            sfd_sample(dimension, sample_size, file_name_csv)
            np.random.seed(iseed) # Set numpy random seed
            latin_hypercube(dimension, sample_size, file_name_csv)
            np.random.seed(iseed) # Set numpy random seed
            sobol_seq(dimension, sample_size, file_name_csv)
            np.random.seed(iseed)
            random_sample(dimension, sample_size, file_name_csv)

def random_sample(dimension, sample_size, file_name_csv):
    best_sample = np.random.random_sample((sample_size, dimension))
    best_dis_random = qmc.discrepancy(best_sample)
    best_maximin_random = maximindist(best_sample)
    best_e_optimality_random = e_optimality(best_sample)
    best_weighted_criteria = np.dot(np.array([best_dis_random, -best_maximin_random, -best_e_optimality_random]), weights)

    for i in range(0,20):
        sample = np.random.random_sample((sample_size, dimension))
        dis_random = qmc.discrepancy(sample)
        maximin_random = maximindist(sample)
        e_optimality_random = e_optimality(sample)
        weighted_criteria = np.dot(np.array([dis_random, -maximin_random, -e_optimality_random]), weights)
        if weighted_criteria < best_weighted_criteria:
            best_sample = sample
            best_weighted_criteria = weighted_criteria
            best_dis_random = dis_random
            best_e_optimality_random = e_optimality_random
            best_maximin_random = maximin_random


    with open(file_name_csv, "a") as file_instance:

        writer = csv.DictWriter(file_instance, fieldnames=fieldnames)
        writer.writerow({
            'method':'Random',
            'discrepancy':best_dis_random,
            'maximin':best_maximin_random,
            'eigen_value':best_e_optimality_random,
            'cumulative':best_weighted_criteria
        })    


def sfd_sample(dimension, sample_size, file_name_csv):
    """ Demonstrate the SamplingModel class on a 2D problem with 10 iterations
        and a desired sample size of 30. """

    import time

    model = sampling_model.SamplingModel( 
        dimension_of_input_space=dimension, 
        grid_size=grid_cells_per_dimension, 
        file_name=file_name_csv,
        no_of_iterations_per_perturbation = no_of_iterations_per_perturbation, # Increase with dimension
        adaptive_sample_size = adaptive_sample_size,
        sample_size = sample_size, # This can be user's choice, e.g., 10
        weights = weights
        )
    start = time.time()
    # Train the model and save results to file_name_csv
    model.initialize()

### Other sampling methods from scipy.qmc ###
              
def latin_hypercube(dimension, sample_size, file_name_csv):
    sampler = qmc.LatinHypercube(dimension)
    sample = sampler.random(sample_size)
            
    dis_lhs = qmc.discrepancy(sample)
    maximin_lhs = maximindist(sample)
    e_optimality_lhs = e_optimality(sample)
    weighted_criteria_lhs = np.dot(np.array([dis_lhs, - maximin_lhs, - e_optimality_lhs]), weights)
    with open(file_name_csv, "a") as file_instance:
        
        writer = csv.DictWriter(file_instance, fieldnames=fieldnames)
        writer.writerow({
            'method':'LHS',
            'discrepancy':dis_lhs,
            'maximin':maximin_lhs,
            'eigen_value':e_optimality_lhs,
            'cumulative':weighted_criteria_lhs
        })       

def sobol_seq(dimension, sample_size, file_name_csv):
    sampler = qmc.Sobol(dimension, scramble=False)
    sample = sampler.random(sample_size)
    dis_sobol = qmc.discrepancy(sample)
    maximin_sobol = maximindist(sample)
    e_optimality_sobol = e_optimality(sample)
    weighted_criteria_sobol = np.dot(np.array([dis_sobol, - maximin_sobol, - e_optimality_sobol]), weights)
    with open(file_name_csv, "a") as file_instance:
        
        writer = csv.DictWriter(file_instance, fieldnames=fieldnames)
        writer.writerow({
            'method':'Sobol',
            'discrepancy':dis_sobol,
            'maximin':maximin_sobol,
            'eigen_value':e_optimality_sobol,
            'cumulative':weighted_criteria_sobol
        })       
        
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

    import sys

    fname = "comparison_dimension_3_seed_"
    for arg in sys.argv[1:]:
        fname = fname + arg
    fname = fname + ".csv"
    for arg in sys.argv[1:]:
        iseed = int(arg)
        comparison(iseed, fname)
