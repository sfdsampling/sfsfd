import csv
from matplotlib import pyplot as plt

# Set templates for Dir and File name
DIR = "../results/"
FNAME = "comparison_dimension_D_seed_S.csv"

# Keys for accessing offsets of methods
method_key = {'SF-SFD': 11,
              'LHS': 12,
              'Sobol': 13,
              'uniform rand': 15}

# Arrays for aggregating results
dims =      [5, 10, 15, 20, 25, 30]
avg_h_sfd = [0, 0,  0,  0,  0,  0]
avg_h_lhs = [0, 0,  0,  0,  0,  0]
avg_h_sob = [0, 0,  0,  0,  0,  0]
avg_h_urs = [0, 0,  0,  0,  0,  0]

# Loop over all dimensions, seeds, and sample sizes
for i, d in enumerate(dims):
    for seed in range(10):
        with open(DIR + FNAME.replace("D", str(d)).replace("S", str(seed)), "r") as fp:
            csv_reader = csv.reader(fp)
            for j, row in enumerate(csv_reader):
                # Aggregate results by offset % 11
                if j > 10:
                    if (j - method_key['SF-SFD']) % 11 == 0:
                        avg_h_sfd[i] += float(row[1])
                    elif (j - method_key['LHS']) % 11 == 0:
                        avg_h_lhs[i] += float(row[1])
                    elif (j - method_key['Sobol']) % 11 == 0:
                        avg_h_sob[i] += float(row[1])
                    elif (j - method_key['uniform rand']) % 11 == 0:
                        avg_h_urs[i] += float(row[1])
    # Average (5 samples size * 10 seeds)
    avg_h_sfd[i] /= 50
    avg_h_lhs[i] /= 50
    avg_h_sob[i] /= 50
    avg_h_urs[i] /= 50

# Set plot style
plt.rc('text', usetex=True)
plt.rc('font', size=13)
# Add labels
plt.ylabel("Avg. discrepancy")
plt.xlabel("Dimension")
# Add lines
plt.plot(dims, avg_h_sfd, '-o', label="SF-SFD")
plt.plot(dims, avg_h_lhs, '-d', label="LHS")
plt.plot(dims, avg_h_sob, '-*', label="Sobol")
plt.plot(dims, avg_h_urs, '-^', label="Unif. rand")
plt.legend(loc="upper left")
plt.tight_layout()
# Show or save figures
#plt.show()
plt.yscale("linear")
plt.savefig("discrep_linear.eps", format="eps")
plt.yscale("log")
plt.savefig("discrep_log.eps", format="eps")
