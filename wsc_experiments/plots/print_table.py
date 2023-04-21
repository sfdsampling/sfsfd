import csv

# Set templates for Dir and File name
DIR = "../results/"
FNAME = "comparison_dimension_D_seed_S.csv"

# Keys for accessing offsets of methods
method_key = {'SF-SFD': 11,
              'LHS': 12,
              'Sobol': 13,
              'uniform rand': 15}

# Arrays for aggregating results
nsample   = [100, 200, 300, 400, 500]
dims      = [5, 10, 15, 20, 25, 30]

# Print table headers
print("""
\\begin{tabular}{|c|l|ccccc|l}
\\cline{1-7}
\\multirow{2}{*}{\\textbf{Dimension}} & \\multicolumn{1}{c|}{\\multirow{2}{*}{\\textbf{$~$ Method $~$}}}
                                      & \\multicolumn{5}{c|}{\\textbf{Sample sizes}}  &  \\\\
\\cline{3-7} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{$~~~~~\\mathbf{100} ~~~~~$}
                                      & \\multicolumn{1}{c|}{$~~~~~\\mathbf{200} ~~~~~$}
                                      & \\multicolumn{1}{c|}{$~~~~~~\\mathbf{300} ~~~~~$}
                                      & \\multicolumn{1}{c|}{$~~~~~\\mathbf{400} ~~~~~$} &
                                      $~~~~~\\mathbf{500} ~~~~~$ &  \\\\
\\cline{1-7}""")

# Loop over all dimensions, seeds, and sample sizes
for d in dims:
    avg_h_sfd = [0,   0,   0,   0,   0]
    avg_h_lhs = [0,   0,   0,   0,   0]
    avg_h_sob = [0,   0,   0,   0,   0]
    avg_h_urs = [0,   0,   0,   0,   0]
    for seed in range(10):
        with open(DIR + FNAME.replace("D", str(d)).replace("S", str(seed)), "r") as fp:
            csv_reader = csv.reader(fp)
            ni = 0
            for j, row in enumerate(csv_reader):
                # Aggregate results by offset % 11
                if j > 10:
                    if (j - method_key['SF-SFD']) % 11 == 0:
                        avg_h_sfd[ni] += float(row[1])
                    elif (j - method_key['LHS']) % 11 == 0:
                        avg_h_lhs[ni] += float(row[1])
                    elif (j - method_key['Sobol']) % 11 == 0:
                        avg_h_sob[ni] += float(row[1])
                    elif (j - method_key['uniform rand']) % 11 == 0:
                        avg_h_urs[ni] += float(row[1])
                        ni += 1
    # Average (10 seeds)
    for i in range(len(nsample)):
        avg_h_sfd[i] /= 10
        avg_h_lhs[i] /= 10
        avg_h_sob[i] /= 10
        avg_h_urs[i] /= 10
    # Print row of table
    print(f"""
\\multirow{{4}}{{*}}{{{d}}}
    & SF-SFD & \\multicolumn{{1}}{{c|}}{{{avg_h_sfd[0]:.4f}}}     & \\multicolumn{{1}}{{c|}}{{{avg_h_sfd[1]:.4f}}}
             & \\multicolumn{{1}}{{c|}}{{{avg_h_sfd[2]:.4f}}}     & \\multicolumn{{1}}{{c|}}{{{avg_h_sfd[3]:.4f}}}
             & {avg_h_sfd[4]:.4f}     &  \\\\ \\cline{{2-7}}
    & LHS    & \\multicolumn{{1}}{{c|}}{{{avg_h_lhs[0]:.4f}}}     & \\multicolumn{{1}}{{c|}}{{{avg_h_lhs[1]:.4f}}}
             & \\multicolumn{{1}}{{c|}}{{{avg_h_lhs[2]:.4f}}}     & \\multicolumn{{1}}{{c|}}{{{avg_h_lhs[3]:.4f}}}
             & {avg_h_lhs[4]:.4f}     &  \\\\ \\cline{{2-7}}
    & Sobol  & \\multicolumn{{1}}{{c|}}{{{avg_h_sob[0]:.4f}}}     & \\multicolumn{{1}}{{c|}}{{{avg_h_sob[1]:.4f}}}
             & \\multicolumn{{1}}{{c|}}{{{avg_h_sob[2]:.4f}}}     & \\multicolumn{{1}}{{c|}}{{{avg_h_sob[3]:.4f}}}
             & {avg_h_sob[4]:.4f}     &  \\\\ \\cline{{2-7}}
    & Unif.\ Rand.\ & \\multicolumn{{1}}{{c|}}{{{avg_h_urs[0]:.4f}}}     & \\multicolumn{{1}}{{c|}}{{{avg_h_urs[1]:.4f}}}
                    & \\multicolumn{{1}}{{c|}}{{{avg_h_urs[2]:.4f}}}     & \\multicolumn{{1}}{{c|}}{{{avg_h_urs[3]:.4f}}}
                    & {avg_h_urs[4]:.4f}     &  \\\\ \\cline{{1-7}}
""")
print("\\end{tabular}")
