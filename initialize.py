import sampling_model2 as samplingModel
import math


dimension = 4
grid_cells_per_dimension = 3
iterations = 100
sample_size = 30
file_name = 'lhc2__dimension_'+str(dimension)+'_' +'grid_cells_'+ str(grid_cells_per_dimension)+'.txt'
model = samplingModel.Model( 
    dimension_of_input_space=dimension, 
    grid_size=grid_cells_per_dimension, 
    file_name=file_name,
    no_of_iterations_per_perturbation = iterations, #Needs to increase with the dimension
    sample_size = sample_size # This can be 10
    )

model.initialize()    




