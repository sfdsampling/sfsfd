import numpy as np
import logging
import scipy.optimize as optimize
from scipy.spatial.distance import pdist
from scipy.stats import qmc
from scipy import fft
from .utils import polar_to_fourier, fourier_to_polar
import csv

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
                 file_name='sfd_sample.txt',
                 no_of_iterations_per_perturbation=100,
                 adaptive_sample_size=0,
                 sample_size=10,
                 bb_budget=None,
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
                we will save the sample history. (defaults to 'sfd_sample.txt')

            no_of_iterations_per_perturbation (int, optional): The number
                of samples that we will take in each iteration of the
                optimization algorithm. (defaults to 100)

            adaptive_sample_size (int, optional): Increase the sample
                size every time this many iterations have elapsed.
                A zero value (default) results in no adaptive increase.

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
        self.sample_size = sample_size
        self.sample_obtained = []
        self.criteria_value_of_sample = (float(self.dimension_of_input_space) *
                                         float(self.dimension_of_input_space) *
                                         np.sqrt(float(self.sample_size)))
        self.discrepancy_value_of_sample = 1.0
        self.min_eigen_value_of_sample = float(self.dimension_of_input_space)
        self.maximin_dist_value_of_sample = float(self.dimension_of_input_space) * np.sqrt(float(self.sample_size))
        self.final_exp_criteria = (float(self.dimension_of_input_space) *
                                   float(self.dimension_of_input_space) *
                                   np.sqrt(float(self.sample_size)))
        self.final_exp_disc = 1.0
        self.final_exp_e_optimality = float(self.dimension_of_input_space)
        self.final_exp_maximin = float(self.dimension_of_input_space) * np.sqrt(float(self.sample_size))
        self.no_of_perturbations_performed = 0.0
        self.no_of_iterations_per_perturbation = \
                                    no_of_iterations_per_perturbation
        self.iterate = 0.0
        self.adaptive_sample_size = adaptive_sample_size
        self.criteria_array_for_all_perturbations = np.array([])
        self.history = []
        if weights is None:
            self.weights = np.ones(3)
        else:
            self.weights = weights
        if bb_budget is None:
            self.bb_budget = 100 * self.grid_size
        else:
            self.bb_budget = bb_budget

    def initialize(self):
        """ Initialize/reset and re-train the sampler using COBYLA.

        Call this method before training a new sampler, or to reset an
        existing sampler.

        """

        # Write status to the output file
        with open(self.file_name, "a") as file_instance:
            file_instance.write("\nThe dimension of the input space is: " +
                                f"{self.dimension_of_input_space}\n")
            
        # Log meta-data
        logging.info("New sample:")
        logging.info(f"  p dimension: {self.dimension_of_input_space}")
        logging.info(f"  sample size: {self.sample_size}")
        logging.info(f"  n grid cell: {self.grid_size}")

        # Tick
        import time
        start = time.time()
        ## Initialize variables
        
        # Start from a uniform distribution
        root_probability_initial_sample = np.array([np.sqrt(1.0 / self.grid_size)
                                           for i in range(self.grid_size)])
        coeff_array = self.fourier_transform(
                                root_prob_mat=root_probability_initial_sample)
        coeff_array_adjusted = []
        for i in coeff_array:
            coeff_array_adjusted.append(np.real(i))
            coeff_array_adjusted.append(np.imag(i))
        angle_array = fourier_to_polar(coeff_array_adjusted)
        # Optimize with COBYLA solver
        optimal_angles_data= optimize.minimize(
            fun = self.iterative_step, 
            x0 = angle_array, 
            method = 'COBYLA',
            options={'maxiter': self.bb_budget}
        )
        # Tock
        walltime = time.time() - start
        # Write history data to the output file
        fieldnames=["method","discrepancy","maximin","eigen_value","cumulative","time_to_sol"]
        with open(self.file_name, "a") as file_instance:
            writer = csv.DictWriter(file_instance, fieldnames=fieldnames)
            writer.writerow({
                'method':'Sf-sfd',
                'discrepancy':self.discrepancy_value_of_sample,
                'maximin':self.maximin_dist_value_of_sample,
                'eigen_value':self.min_eigen_value_of_sample,
                'cumulative':self.criteria_value_of_sample,
                'time_to_sol':walltime
            })

    def generate_initial_sample(self):
        """ Create an initial sample.

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
            initial_sample (numpy.ndarray): An initial numpy array of sample 
                points of size dim X num coeffs that is used to generate our 
                initial Fourier coeffs.

        Returns:
            numpy.ndarray: Flattened grid corresponding  

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
            initial_sample_discretized (numpy.ndarray):binned samples

        Returns:
            numpy.ndarray: pdf for the binned array

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
            root_prob_mat (numpy.ndarray): A numpy array of equally spaced 
                finite sample points

        Returns:
            numpy.ndarray: Fourier coefficients corresponding to DFT of
                root_prob_mat

        """

        cm_root = fft.fft(root_prob_mat)
        numpy_cm_root = np.array(cm_root) / np.sqrt(self.grid_size)
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
            fourier_root_cm.append(complex(fourier_root_cm_flat[2*i],
                                           fourier_root_cm_flat[(2*i)+1]))
        fourier_root_cm = np.asarray(fourier_root_cm)
        probability_matrix_obtained = self.create_prob_distribution_1D(
                                                            fourier_root_cm)
        # sum of all entries of prob matrix is 1
        criteria_array_for_iterations_in_perturbation = []
        discrepancy_array_for_iterations_in_perturbation = []
        maximin_dist_array_for_iterations_in_perturbation = []
        e_optimality_array_for_iterations_in_perturbation = []
        optimal_sample = []
        optimal_sample_criteria_value = np.infty

        for iteration in range(self.no_of_iterations_per_perturbation):
            sample_created = self.sampling_from_iid_distribution(
                                                probability_matrix_obtained)
            criteria_value_data = self.criteria_result(sample_created)
            criteria_value = criteria_value_data[0]
            
            if criteria_value <= optimal_sample_criteria_value:
                optimal_sample = sample_created
                optimal_sample_criteria_value = criteria_value
                optimal_sample_disc_value = criteria_value_data[1]
                optimal_sample_maximin_value = criteria_value_data[2]
                optimal_sample_eigen_value = criteria_value_data[3]
            
            criteria_array_for_iterations_in_perturbation.append(
                                                criteria_value)
            maximin_dist_array_for_iterations_in_perturbation.append(
                                                criteria_value_data[2])
            discrepancy_array_for_iterations_in_perturbation.append(
                                                criteria_value_data[1])
            e_optimality_array_for_iterations_in_perturbation.append(
                                                criteria_value_data[3])
        # Update the sample size according to adaptive schedule
        self.iterate += 1
        if self.adaptive_sample_size:
            if self.iterate % self.adaptive_sample_size == 0:
                self.no_of_iterations_per_perturbation += 1

        # Calcualte expected values
        criteria_value_for_the_perturbation = \
                    (sum(criteria_array_for_iterations_in_perturbation) /
                     len(criteria_array_for_iterations_in_perturbation))
        disc_value_for_the_perturbation = \
                (sum(discrepancy_array_for_iterations_in_perturbation) /
                len(discrepancy_array_for_iterations_in_perturbation))
        maximin_value_for_the_perturbation = \
                (sum(maximin_dist_array_for_iterations_in_perturbation) /
                len(maximin_dist_array_for_iterations_in_perturbation))
        eigenvalue_value_for_the_perturbation = \
                (sum(e_optimality_array_for_iterations_in_perturbation) /
                len(e_optimality_array_for_iterations_in_perturbation))
        
        
        self.criteria_array_for_all_perturbations = np.append(
                                self.criteria_array_for_all_perturbations,
                                criteria_value_for_the_perturbation)

        # If this is an improvement, record best and average-case
        if criteria_value_for_the_perturbation < self.final_exp_criteria:
            self.sample_obtained = optimal_sample
            self.criteria_value_of_sample = optimal_sample_criteria_value
            self.discrepancy_value_of_sample = optimal_sample_disc_value
            self.min_eigen_value_of_sample = optimal_sample_eigen_value
            self.maximin_dist_value_of_sample = optimal_sample_maximin_value

            self.final_exp_criteria = criteria_value_for_the_perturbation
            self.final_exp_disc = disc_value_for_the_perturbation
            self.final_exp_maximin = maximin_value_for_the_perturbation
            self.final_exp_e_optimality = eigenvalue_value_for_the_perturbation

        logging.info(f"d: {self.dimension_of_input_space}, " +
                     f"n: {self.sample_size}, " +
                     f"iter: {self.iterate}")
        logging.info(f"  sample size: {self.no_of_iterations_per_perturbation}")
        logging.info(f"  discrep: {optimal_sample_disc_value}")
        logging.info(f"  e-optim: {optimal_sample_eigen_value}")
        logging.info(f"  maximin: {optimal_sample_maximin_value}")
        logging.info(f"  normalized score: {criteria_value_for_the_perturbation}")

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

        if self.weights[0] > 0:
            discrepancy = qmc.discrepancy(sample)
        else:
            discrepancy = 1.0
        if self.weights[1] > 0:
            maximindistance = np.min(pdist(sample)) # By default Euclidean distance
        else:
            maximindistance = 0.0
        if self.weights[2] > 0:
            s = np.linalg.svd(sample, compute_uv=False)
            min_eigenvalue = np.min(s)
        else:
            min_eigenvalue = 0.0
        # Record all 3 samples in the history array
        self.history.append([discrepancy, maximindistance, min_eigenvalue])
        # Calculate weighted average
        b1 = np.sqrt(float(self.dimension_of_input_space))
        b2 = float(self.dimension_of_input_space) * np.sqrt(float(self.sample_size))
        result = np.prod(np.array([discrepancy**self.weights[0],
                                   (b1 - maximindistance) ** self.weights[1],
                                   (b2 - min_eigenvalue) ** self.weights[2]]))
        return [result, discrepancy, maximindistance, min_eigenvalue]

    def create_prob_distribution_1D(self, fourier_root_cm):
        ''' Generate a 1D PDF from a 1D list of sqrt-Fourier coeffs.

        Args:
            fourier_root_cm (numpy.ndarray): The sqrt-Fourier coeffs for a 1D
                PDF, from which we calculate the dD PDF values assuming iid.

        Returns:
            numpy.ndarray: PDF values caluclated

        '''

        prob_ifft_squareroot = fft.ifft(fourier_root_cm * np.sqrt(self.grid_size))
        
        ## We can assert this
        # Sampling to find the expected value of discrepancy
        # Find the probability matrix by squaring and taking absolute
        #         value of squareroot FT we obtained
        # For numerical reasons, take the real part and re-normalize by sum
        prob_matrix_1d = pow(np.abs(prob_ifft_squareroot),2)
        prob_matrix_1d = prob_matrix_1d / np.sum(prob_matrix_1d)
        return prob_matrix_1d

    def create_prob_distribution(self, fourier_root_cm):
        ''' Generate a PDF from a list of sqrt-Fourier coeffs.

        Args:
            fourier_root_cm (numpy.ndarray): The sqrt-Fourier coeffs for which
                we calculate the PDF values.

        Returns:
            numpy.ndarray: PDF values caluclated

        '''

        prob_ifft_squareroot = np.array(fft.ifft(fourier_root_cm))
        ## We can assert this
        #file_instance.write(f"The sum of probabilities is: {sum_l2_prob}\n")
        # Sampling to find the expected value of discrepancy
        # Find the probability matrix by squaring and taking absolute
        #         value of squareroot FT we obtained
        probability_matrix_obtained = abs(pow(prob_ifft_squareroot,2))
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
        samples_cell_mapping = np.array([])
        # Map the samples to the Random variable
        for each_sample in our_sample:
            diff = each_sample
            i=0
            while diff >0:
                diff -= probability_distribution[i]
                i+=1
            samples_cell_mapping = np.append(samples_cell_mapping,i-1)
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

    def sampling_from_iid_distribution(self, probability_distribution):
        """ Generate a sample from a PDF.

        Given a PDF, generate a random sample of self.sample_size points.

        Args:
            probability_distribution (numpy.ndarray): A 1D distribution
                representation of the PDF values.

        Returns:
            numpy.ndarray: A dD sample from probability_distribution, assuming
            all d dimensions are iid, of size self.sample_size.

        """

        our_sample = np.random.rand(self.sample_size *
                                    self.dimension_of_input_space)
        samples_cell_mapping = np.array([])
        # Map the samples to the Random variable
        for each_sample in our_sample:
            diff = each_sample
            i=0
            while diff >0:
                diff -= probability_distribution[i]
                i+=1
            samples_cell_mapping = np.append(samples_cell_mapping,np.array([i-1]))
        samples_cell_mapping = samples_cell_mapping.reshape((self.sample_size,
                                                  self.dimension_of_input_space))
        sample_created = None
        # These coordinates are (x_d, ..., x_1)
        for each_sample_point in samples_cell_mapping:
            d = self.dimension_of_input_space
            high=[]
            low=[]
            for i in range(self.dimension_of_input_space):
                q = each_sample_point[i]
                low.append(q*self.delta)
                high.append((q+1)*self.delta)
            sample_point = np.random.uniform(low=low, high=high,
                                        size=(1,self.dimension_of_input_space))
            if sample_created is None:
                sample_created = sample_point.copy()
            else:
                sample_created = np.append(sample_created, sample_point, axis=0)
        return sample_created
