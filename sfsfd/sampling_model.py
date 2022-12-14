import numpy as np
import scipy.optimize as optimize
from scipy.spatial.distance import pdist
from scipy.stats import qmc
from scipy import fft
from .utils import polar_to_fourier, fourier_to_polar


class SamplingModel:
    """ This class generates sampling models via optimization of Fourier coefs.

    The public methods are:
     - __init__ : This is the constructor
     - initialize : This resets and trains a sampler object with COBYLA

    The private methods are:
     - generate_initial_sample : This generates an initial sample to tune
     - discretization_of_points : This "bins" the samples
     - grid_to_cell_mapping_probability : This maps the binned sample to a PDF
     - fourier_transform : Computes the Fourier transform of the sqrt p-vals
     - iterative_step : Generate an objective score for a given input array
     - criteria_result : Calculate (weighted) average of all 3 criteria
     - create_prob_distribution : Given a set of Fourier roots, calculate PDF
     - sampling_from_distribution : Given a PDF, generate a sample

    """

    def __init__(self, 
                 dimension_of_input_space=2, 
                 grid_size=3, 
                 file_name='sample.txt',
                 no_of_iterations_per_perturbation=100,
                 sample_size=10,
                 weights=None):
        """ Constructor for the SamplingModel class.

        Args:
            dimension_of_input_space (int, optional): The dimension of the
                space to sample. (defaults to 2 when omitted)

            grid_size (int, optional): The number of cells per dimension
                in the discretized PDF. Note: the dimension of the
                optimization problem grows like O(grid_size ^ d).
                (defaults to 3)

            file_name (str, optional): The name of the output file where.
                we will save the sample history. (defaults to 'sample.txt')

            no_of_iterations_per_perturbation (int, optional): The number
                of samples that we will take in each iteration of the
                optimization algorithm. (defaults to 100)

            sample_size (int, optional): The size of the sample that we
                would like to draw from our trained sampler. (defaults to 10)

            weights (numpy.ndarray, optional): How to weight the 3 criteria
                -- L2 discrepancy, maximin distance, and E-optimality --
                when calculating objective scores. Defaults to equal
                weighting np.array([0.33, 0.33, 0.33])

        """

        self.dimension_of_input_space = dimension_of_input_space
        self.grid_size = grid_size
        self.delta = 1/self.grid_size
        self.no_of_coefficients = pow(self.grid_size,
                                      self.dimension_of_input_space)
        self.file_name = file_name
        self.coeff_array = []
        self.probability_matrix =[]
        self.optimal_sample = np.array([])
        self.disc_of_optimal_sample = 1
        self.sample_size = sample_size
        self.sample_obtained = []
        self.criteria_value_of_sample = 1
        self.no_of_perturbations_performed = 0
        self.no_of_iterations_per_perturbation = \
                                    no_of_iterations_per_perturbation
        self.criteria_array_for_all_perturbations = np.array([])
        self.history = []
        if weights is None:
            self.weights = np.ones(3) / 3.0
        else:
            self.weights = weights

    def initialize(self):
        """ Initialize/reset and re-train the sampler using COBYLA.

        Call this method before training a new sampler, or to reset an
        existing sampler.

        """

        # Write status to the output file
        with open(self.file_name, "a") as file_instance:
            file_instance.write("The dimension of the input space is: " +
                                f"{self.dimension_of_input_space}\n")
            file_instance.write("Number of grid cells per dimension = " +
                                f"{self.grid_size}\n")
            file_instance.write("Total number of fourier coefficients = " +
                                f"{self.no_of_coefficients}\n")
        # Initialize variables
        initial_sample = self.generate_initial_sample()
        discretized_sample = self.discretization_of_points(
                                initial_sample=initial_sample)
        probability_initial_sample = self.grid_to_cell_mapping_probability(
                                initial_sample_discretized=discretized_sample)
        root_probability_initial_sample = np.sqrt(probability_initial_sample)
        coeff_array = self.fourier_transform(
                                root_prob_mat=root_probability_initial_sample)
        # Create an initial angle array (x0) for optimization with COBYLA
        coeff_array_adjusted = []
        for i in coeff_array:
            coeff_array_adjusted.append(np.real(i))
            coeff_array_adjusted.append(np.imag(i))
        angle_array = fourier_to_polar(coeff_array_adjusted)
        # Optimize with COBYLA solver
        optimal_angles_data= optimize.minimize(
            fun = self.iterative_step, 
            x0 = angle_array, 
            method = 'COBYLA'
        )
        # Write history data to the output file
        with open(self.file_name, "a") as file_instance:
            file_instance.write("Number of perturbations = " +
                                f"{self.no_of_perturbations_performed}\n")
            file_instance.write("The sample created is: " +
                                f"{self.sample_obtained}\n")
            file_instance.write("Criteria value of the sample is: " +
                                f"{self.criteria_value_of_sample}\n")
            file_instance.write("Expected value of discrepancy of the sample:"
                                + f" {optimal_angles_data['fun']}\n")
            file_instance.write("Criteria value for all perturbations: " +
                        f"{repr(self.criteria_array_for_all_perturbations)}\n")
        
    def generate_initial_sample(self):
        """ Step 1: Create an initial sample.

        Returns:

            numpy.ndarray : An initial numpy array of sample points of size
                dim X num coeffs that can be used to generate our initial
                Fourier coeffs.

        """

        initial_sample = qmc.LatinHypercube(
                    d=self.dimension_of_input_space).random(
                                                    self.no_of_coefficients)
        ## Uncomment to dump additional diagnostics to output file
        #with open(self.file_name, "a") as file_instance:
        #    file_instance.write("Initial sample by Latin Hypercube " +
        #                        f"Sequence is: {initial_sample}\n")
        #    file_instance.write("This sample is used to create the " +
        #                        "Fourier series representation of the " +
        #                        "distribution function\n")
        #    file_instance.write("Length of the Fourier sample: " +
        #                        f"{len(initial_sample)}\n")
        return initial_sample

    def discretization_of_points(self, initial_sample):
        ''' Generate a discretization of the input sample.

        This also handles samples with dimension 1.
        We consider their coordinate as grid_size-1th coordinate

        Args:
            initial_sample (numpy.ndarray): TODO

        Returns:
            numpy.ndarray: TODO

        '''

        initial_sample_discretized = []
        for each_point in initial_sample:
            disc_point = [(np.floor(each_point[i]/self.delta) - 1)
                          if (np.floor(each_point[i]/self.delta) ==
                              self.grid_size)
                          else (np.floor(each_point[i]/self.delta))
                          for i in range(self.dimension_of_input_space)]
            initial_sample_discretized.append(disc_point)
        return initial_sample_discretized

    def grid_to_cell_mapping_probability(self, initial_sample_discretized):
        """ Convert binned samples from above into PDF values.

        Args:
            initial_sample_discretized (numpy.ndarray): TODO

        Returns:
            numpy.ndarray: TODO

        """

        pmf = np.zeros(self.no_of_coefficients)
        for each_sample in initial_sample_discretized:
            pmf_index = 0
            for i in range(self.dimension_of_input_space):
                pmf_index += int(each_sample[i] * pow(self.grid_size, i))
            pmf[pmf_index] += 1
        return pmf / len(initial_sample_discretized)
    
    def fourier_transform(self, root_prob_mat):
        """ Convert root-PDF values into Fourier coefficients of sqrt(PDF).

        Args:
            root_prob_mat (numpy.ndarray): TODO

        Returns:
            numpy.ndarray: TODO

        """

        cm_root = fft.fft(root_prob_mat)
        numpy_cm_root = np.array(cm_root)
        numpy_cm_root = numpy_cm_root/np.sqrt(self.no_of_coefficients)
        return numpy_cm_root

    def iterative_step(self, angle_array):
        """ Calculate a (weighted) criteria score for the given angle array.

        Args:
            angle_array (numpy.ndarray): Polar coordinate transform of complex
                Fourier coeffs on the Bloch sphere.

        Returns:
            float: The (weighted) average of our 3 performance criteria

        """

        # Incrememnt our iteration counter
        self.no_of_perturbations_performed += 1
        # Convert to PDF
        fourier_root_cm_flat = polar_to_fourier(angle_array)
        index = int((len(fourier_root_cm_flat))/2)
        fourier_root_cm = []
        for i in range(index):
            fourier_root_cm.append([complex(fourier_root_cm_flat[2*i],
                                            fourier_root_cm_flat[(2*i)+1])])
        
        probability_matrix_obtained = self.create_prob_distribution(
                                                            fourier_root_cm)
        # sum of all entries of prob matrix is 1
        criteria_array_for_iterations_in_perturbation = []
        optimal_sample = []
        optimal_sample_criteria_value = 1
        for iteration in range(self.no_of_iterations_per_perturbation):
            sample_created = self.sampling_from_distribution(
                                                probability_matrix_obtained)
            criteria_value = self.criteria_result(sample_created)
            if criteria_value <= optimal_sample_criteria_value:
                optimal_sample = sample_created
                optimal_sample_criteria_value = criteria_value
            criteria_array_for_iterations_in_perturbation.append(
                                                                criteria_value)
        criteria_value_for_the_perturbation = \
                    (sum(criteria_array_for_iterations_in_perturbation) /
                     len(criteria_array_for_iterations_in_perturbation))
        self.criteria_array_for_all_perturbations = np.append(
                                self.criteria_array_for_all_perturbations,
                                criteria_value_for_the_perturbation)
        self.sample_obtained = optimal_sample
        self.criteria_value_of_sample = optimal_sample_criteria_value
        return criteria_value_for_the_perturbation

    def criteria_result(self, sample):
        ''' Generate a criteria score that can be minimized for the given
        sample, based on following 3 criteria:
         1. We want to maximize the maximin distance
         2. We want as low discrepancy as possible
         3. We want E-optimality criteria

        Args:
            sample (numpy.ndarray): The sample for which we calculate the
                criteria scores.

        Returns:
            float: The (weighted) average of 3 criteria above

        '''

        maximindistance = max(pdist(sample)) # By default Euclidean distance
        discrepancy = qmc.discrepancy(sample)
        sample_arr = np.array([arr.tolist() for arr in sample])
        t = sample_arr.T # 4x10
        u,s,v = np.linalg.svd(t)
        min_eigenvalue = np.min(s)
        # Record all 3 samples in the history array
        self.history.append([discrepancy, maximindistance, min_eigenvalue])
        # Calculate weighted average
        result = np.dot(self.weights,
                        np.array([discrepancy, -maximindistance,
                                  -min_eigenvalue]))
        #if result < self.disc_of_optimal_sample:
        #    self.optimal_sample = sample
        #    self.disc_of_optimal_sample = result
        return result

    def create_prob_distribution(self, fourier_root_cm):
        ''' Generate a PDF from a list of sqrt-Fourier coeffs.

        Args:
            fourier_root_cm (numpy.ndarray): The sqrt-Fourier coeffs for which
                we calculate the PDF values.

        Returns:
            numpy.ndarray: PDF values caluclated

        '''

        prob_ifft_squareroot = np.array(fft.ifft(fourier_root_cm))
        '''
        sum_l2_prob = 0
        for i in prob_ifft_squareroot:
            sum_l2_prob+= pow(abs(i), 2)
        '''
        ## We can assert this
        #file_instance.write(f"The sum of probabilities is: {sum_l2_prob}\n")
        # Sampling to find the expected value of discrepancy
        # There are total 6 steps
        # Step 1: Find the probability matrix by squaring and taking absolute
        #         value of squareroot FT we obtained
        probability_matrix_obtained = abs(pow(prob_ifft_squareroot,2))
        #file_instance.write("The probability matrix obtained is: " +
        #                    "probability_matrix_obtained\n")
        return probability_matrix_obtained

    def sampling_from_distribution(self, probability_distribution):
        """ Generate a sample from a PDF.

        Given a PDF, generate a random sample of self.sample_size points.

        Args:
            probability_distribution (numpy.ndarray): A representation of the
                PDF values.

        Returns:
            numpy.ndarray: A sample from probability_distribution of size
            self.sample_size.

        """

        our_sample = np.random.rand(self.sample_size)
        #with open(self.file_name, "a") as file_instance:
        #    file_instance.write("The sample to be mapped according to our " +
        #                        f"distribution is: {our_sample}\n")
        samples_cell_mapping = np.array([])
        # Step 3: Map the samples to the Random variable
        for each_sample in our_sample:
            diff = each_sample
            i=0
            while diff >0:
                diff -= probability_distribution[i]
                i+=1
            samples_cell_mapping = np.append(samples_cell_mapping,i-1)
        #file_instance.write(samples_cell_mapping)
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
            assert(len(low) == len(grid_coordinates) ==
                   self.dimension_of_input_space)
            grid_map_sample.append(grid_coordinates)
            sample_point = np.random.uniform(low=low, high=high,
                                        size=(1,self.dimension_of_input_space))
            sample_created.append(sample_point[0])
        return sample_created
