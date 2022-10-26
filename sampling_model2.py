#from matplotlib.pyplot import grid
from unittest import result
import numpy as np
from scipy.stats import qmc
import math
from scipy import fft
from numpy import imag, random, real
import scipy.optimize as optimize


def polar_to_fourier(theta):
    """ Convert polar angles into Fourier (wave function) coefficients.

    Args:
        theta (np.ndarray): A numpy array of length n - 1 containing polar
            angles.
    Returns:
        np.ndarray: A numpy array of length n containing the Fourier
        coefficients corresponding to the given polar angles.
    """
    # Initialize output array
    n = len(theta) + 1
    cx = np.ones(n)
    # Calculate the output array
    for i in range(n):
        for j in range(n - 1 - i):
            cx[i] *= np.cos(theta[j])
        if i > 0:
            cx[i] *= np.sin(theta[n - 1 - i])
    return cx

def fourier_to_polar(cx):
    """ Convert Fourier (wave function) coefficients into polar angles.
    Args:
        cx (np.ndarray): A numpy array of length n containing the Fourier
            coefficients corresponding to the given polar angles.

    Returns:
        theta (np.ndarray): A numpy array of length n - 1 containing polar
        angles.

    """
    # Initialize output array
    n = len(cx)
    theta = np.ones(n - 1)
    # Calculate the output array
    for i in range(n - 1):
        xx = cx[n - 1 - i]
        for j in range(i):
            xx /= np.cos(theta[j])
        theta[i] = np.arcsin(xx)
    return theta

if __name__ == "__main__":
    # Embed and extract a diagonal vector
    x1 = np.ones(6) / np.sqrt(6.0)
    theta1 = fourier_to_polar(x1)
    assert(np.linalg.norm(polar_to_fourier(theta1) - x1) < 1.0e-8)
    # Extract and embed the first basis vector
    theta2 = np.ones(3) * 2.0 * np.pi
    x2 = polar_to_fourier(theta2)
    assert(np.linalg.norm(x2 - np.eye(4)[0]) < 1.0e-8)
    assert(np.linalg.norm(theta2 - fourier_to_polar(x2)) < 1.0e-8 or
           np.linalg.norm(fourier_to_polar(x2)) < 1.0e-8)
    # Extract and embed the first basis vector again
    theta3 = np.zeros(11)
    x3 = polar_to_fourier(theta3)
    assert(np.linalg.norm(x3 - np.eye(12)[0]) < 1.0e-8)
    assert(np.linalg.norm(theta3 - fourier_to_polar(x3)) < 1.0e-8 or
           np.linalg.norm(np.ones(11) * 2.0 * np.pi - fourier_to_polar(x2))
           < 1.0e-8)
# No of perturbation for 2-dimensional case = 

class Model:
    def __init__(
        self, 
        dimension_of_input_space=2, 
        grid_size=3, 
        file_name='sample.txt',
        no_of_iterations_per_perturbation = 100,
        sample_size = 10
    ):
        
        self.dimension_of_input_space = dimension_of_input_space
        self.grid_size = grid_size
        self.delta = 1/self.grid_size
        self.no_of_coefficients = pow(self.grid_size, self.dimension_of_input_space)
        self.file_name = file_name
        self.coeff_array = []
        self.probability_matrix =[]
        self.optimal_sample = np.array([])
        self.disc_of_optimal_sample = 1
        self.sample_size = sample_size
        self.sample_obtained = []
        self.criteria_value_of_sample = 1
        self.no_of_perturbations_performed = 0
        self.no_of_iterations_per_perturbation = no_of_iterations_per_perturbation
        #self.discrepancy_arr_for_all_perturbation = []
        self.criteria_array_for_all_perturbations = np.array([])
        pass
    
    def initialize(self):
    # initial calculations

        file_instance = open(self.file_name, "a")
        print("The dimension of the input space is: ", self.dimension_of_input_space, file=file_instance)
        print("Number of grid cells per dimension = ", self.grid_size, file=file_instance)
        print("Total number of fourier coefficients = ", self.no_of_coefficients, file=file_instance)
        file_instance.close()

        initial_sample = self.generate_initial_sample()
        discretized_sample = self.discretization_of_points(initial_sample=initial_sample)
        probability_initial_sample = self.grid_to_cell_mapping_probability(initial_sample_discretized=discretized_sample)
        root_probability_initial_sample = np.sqrt(probability_initial_sample)
        coeff_array = self.fourier_transform(root_prob_mat=root_probability_initial_sample)
        
        
        # Create an angle array for optimization with COBYLA
        coeff_array_adjusted = []
        for i in coeff_array:
            coeff_array_adjusted.append(real(i))
            coeff_array_adjusted.append(imag(i))
        angle_array = fourier_to_polar(coeff_array_adjusted)
        #self.iterative_step(angle_array)

        optimal_angles_data= optimize.minimize(
            fun = self.iterative_step, 
            x0 = angle_array, 
            method = 'COBYLA'
        )
        
        file_instance = open(self.file_name, "a")
        print("Number of perturbations = ", self.no_of_perturbations_performed, file=file_instance)
        print("The sample created is: ", self.sample_obtained, file=file_instance )
        print("Discrepancy of the sample is: ", self.criteria_value_of_sample, file=file_instance)
        print("Expected value of discrepancy of the sample: ", optimal_angles_data['fun'], file=file_instance)
        print("Criteria value for all perturbations", repr(self.criteria_array_for_all_perturbations), file=file_instance)
        file_instance.close()
        
    # Step 1: Create an initial sample
    def generate_initial_sample(self):
        initial_sample =  qmc.LatinHypercube(d =self.dimension_of_input_space).random(self.no_of_coefficients)
        
        file_instance = open(self.file_name, "a")
        print("Initial sample by Latin Hypercube Sequence is:", initial_sample, file=file_instance)
        print("This sample is used to create the fourier series representation of the distribution function", file=file_instance)
        print("Length of the sample: ", len(initial_sample),file=file_instance)
        file_instance.close()
        
        return initial_sample
   
    # Checked, this is correct
    def discretization_of_points(self, initial_sample):
        '''
        This also takes care for samples with coordinates 1.
        We consider their coordinate as grid_size-1th coordinate
        '''
        report_file = open("report_file_data.txt", "a")
        print("initial sample by halton sequence:", initial_sample, file=report_file)
        initial_sample_discretized = []
        for each_point in initial_sample:
            disc_point = [(math.floor(each_point[i]/self.delta) - 1) if math.floor(each_point[i]/self.delta) == self.grid_size else (math.floor(each_point[i]/self.delta)) for i in range(self.dimension_of_input_space)]
            initial_sample_discretized.append(disc_point)
        
        print("Discretization of the sample", initial_sample_discretized, file=report_file)
        report_file.close()
        return initial_sample_discretized

    def grid_to_cell_mapping_probability(self, initial_sample_discretized):
        report_file = open("report_file_data.txt", "a")
        pmf=np.zeros(self.no_of_coefficients)
        for each_sample in initial_sample_discretized:
            pmf_index = 0
            for i in range(self.dimension_of_input_space):
                pmf_index += each_sample[i] * pow(self.grid_size, i)
            pmf[pmf_index] +=1
        print("pmf is: ",pmf, file=report_file)
        report_file.close()
        return pmf/len(initial_sample_discretized)
    
    def fourier_transform(self, root_prob_mat):
        cm_root = fft.fft(root_prob_mat)
        numpy_cm_root = np.array(cm_root)
        numpy_cm_root = numpy_cm_root/math.sqrt(self.no_of_coefficients)
        return numpy_cm_root

    def iterative_step(self, angle_array):
        
        self.no_of_perturbations_performed += 1
        
        fourier_root_cm_flat = polar_to_fourier(angle_array)
        index = int((len(fourier_root_cm_flat))/2)
        fourier_root_cm = []
        for i in range(index):
            fourier_root_cm.append([complex(fourier_root_cm_flat[2*i],fourier_root_cm_flat[(2*i)+1])])
        
        probability_matrix_obtained = self.create_prob_distribution(fourier_root_cm)
        # sum of all entries of prob matrix is 1
        criteria_array_for_iterations_in_perturbation = []
        optimal_sample = []
        optimal_sample_criteria_value = 1
        for iteration in range(self.no_of_iterations_per_perturbation):
            sample_created = self.sampling_from_distribution(probability_matrix_obtained)
            criteria_value = self.criteria_result(sample_created)
            if criteria_value <= optimal_sample_criteria_value:
                optimal_sample = sample_created
                optimal_sample_criteria_value = criteria_value
            criteria_array_for_iterations_in_perturbation.append(criteria_value)

        criteria_value_for_the_perturbation = sum(criteria_array_for_iterations_in_perturbation)/len(criteria_array_for_iterations_in_perturbation)
        self.criteria_array_for_all_perturbations = np.append(self.criteria_array_for_all_perturbations,criteria_value_for_the_perturbation)
        self.sample_obtained = optimal_sample
        self.criteria_value_of_sample = optimal_sample_criteria_value
        return criteria_value_for_the_perturbation


    def criteria_result(self, sample):
        result = qmc.discrepancy(sample)
        '''
        if result < self.disc_of_optimal_sample:
            self.optimal_sample = sample
            self.disc_of_optimal_sample = result
        '''
        return result

    def create_prob_distribution(self, fourier_root_cm):

        prob_ifft_squareroot = np.array(fft.ifft(fourier_root_cm))
        '''
        sum_l2_prob = 0
        for i in prob_ifft_squareroot:
            sum_l2_prob+= pow(abs(i), 2)
        '''
        #We can assert this
        #print("The sum of probabilities is: ", sum_l2_prob, file=file_instance)
        # Sampling to find the expected value of discrepancy
        # There are total 6 steps
        # Step 1: Find the probability matrix by squaring and taking absolute value of squareroot FT we obtained
        probability_matrix_obtained = abs(pow(prob_ifft_squareroot,2))
        #print("The probability matrix obtained is: ", probability_matrix_obtained, file=file_instance )
        return probability_matrix_obtained

    def sampling_from_distribution(self, probability_distribution):
        our_sample = random.rand(self.sample_size)
        #file_instance = open(self.file_name, "a")
        #print("The sample to be mapped according to our distribution is: ", our_sample, file=file_instance)
        samples_cell_mapping = np.array([])
            # Step 3: Map the samples to the Random variable
        for each_sample in our_sample:
            diff = each_sample
            i=0
            while diff >0:
                diff -= probability_distribution[i]
                i+=1
            samples_cell_mapping = np.append(samples_cell_mapping,i-1)
        #print(samples_cell_mapping, file=file_instance)    
        
        grid_map_sample = []
        sample_created = []
        # These coordinates are (x_d, ..., x_1)
        for each_sample_point in samples_cell_mapping:
            d = self.dimension_of_input_space
            grid_coordinates= []
            high=[]
            low=[]
            r = each_sample_point
            for i in range(self.dimension_of_input_space-1):
                
                q = r // pow(self.grid_size, d-i-1)
                r = r % pow(self.grid_size, d-i-1)
                grid_coordinates.append(q)
                low.append(q*self.delta)
                high.append((q+1)*self.delta)

            grid_coordinates.append(r)
            high.append((r+1)*self.delta)
            low.append(r*self.delta)
            assert(len(low) ==len(high))
            assert(len(low) == len(grid_coordinates) == self.dimension_of_input_space)
            grid_map_sample.append(grid_coordinates)
            sample_point = np.random.uniform(low=low, high=high, size=(1,self.dimension_of_input_space))
            sample_created.append(sample_point[0])
        return sample_created