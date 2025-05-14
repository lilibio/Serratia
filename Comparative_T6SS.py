import math
import numpy as np
from scipy.stats import chi2

# sed -i '1s/^[^\t]*\t//' output_matrix.txt
# sed -i '1s/^/\t/' output_matrix.txt
# head -n 1 output_matrix.txt > header
# cat output_matrix.txt|grep "*" > main.txt
# cat header main.txt > orthologue_file.txt
# Written by Qingtian Tim Guan, Version 2.1, chi-square value added in this version. Debug the core gene issue which will cause the expected value of a gene absence in t6ss+/- to be 0, which cannot be divided.
# Organize the "orthologue_file.txt", "presence_genomes.txt","absence_genomes.txt" in the same folder and run python3 Comparative_T6SS.py > LogPA.txt
# Read orthologue presence/absence file
with open("orthologue_file.txt", "r") as ortho_file:
    ortho_lines = ortho_file.readlines()

# Extract genome names from the first line of the orthologue file
    #print(ortho_lines[0])
genome_names = ortho_lines[0].strip().split('\t')[1:]

# Read T6SS presence genomes list
with open("presence_genomes.txt", "r") as presence_file:
    t6ss_presence_genomes = presence_file.read().splitlines()

# Read T6SS absence genomes list
with open("absence_genomes.txt", "r") as absence_file:
    t6ss_absence_genomes = absence_file.read().splitlines()

# Initialize lists to store presence and absence counts for each group
t6ss_presence_counts = []
t6ss_absence_counts = []

# Initialize lists to store p-values for each group
p_values = []

# Iterate through each group in the orthologue file
for line in ortho_lines[1:]:
    group, *genes = line.strip().split('\t')
    t6ss_presence_count = sum(1 for genome in t6ss_presence_genomes if any(gene != '*' for gene in genes[genome_names.index(genome)].split(', ')))
    t6ss_absence_count = sum(1 for genome in t6ss_absence_genomes if any(gene != '*' for gene in genes[genome_names.index(genome)].split(', ')))
    t6ss_presence_counts.append(t6ss_presence_count)
    t6ss_absence_counts.append(t6ss_absence_count)

    # Create a 2x2 contingency table for the chi-square test
    observed = np.array([[t6ss_presence_count, t6ss_absence_count],
                         [len(t6ss_presence_genomes) - t6ss_presence_count, len(t6ss_absence_genomes) - t6ss_absence_count]])

    # Calculate the expected frequencies for each cell
    row_totals = np.sum(observed, axis=1)
    col_totals = np.sum(observed, axis=0)
    total = np.sum(observed)
    expected = np.outer(row_totals, col_totals) / total

    # Calculate the chi-square statistic
    chi2_statistic = np.sum((observed - expected)**2 / expected)

    # Calculate the degrees of freedom
    degrees_of_freedom = (observed.shape[0] - 1) * (observed.shape[1] - 1)

    # Calculate the chi-square p-value using the survival function (1 - CDF)
    p = 1 - chi2.cdf(chi2_statistic, degrees_of_freedom)
    p_values.append(p)

# Calculate the values and sort the data based on the value in descending order
sorted_group_data = []
for i, (presence_count, absence_count, p) in enumerate(zip(t6ss_presence_counts, t6ss_absence_counts, p_values), start=1):
    value = math.log2((presence_count * len(t6ss_absence_genomes) / len(t6ss_presence_genomes) + 1) / (absence_count + 1))
    sorted_group_data.append((i, presence_count, absence_count, value, p))

# Sort the data based on the value in descending order
sorted_group_data.sort(key=lambda x: x[3], reverse=True)

# Print the results
for rank, (i, presence_count, absence_count, value, p) in enumerate(sorted_group_data, start=1):
    group_line = ortho_lines[i].strip()
    print(f"Group {i}:\tNo. T6SS+:{presence_count}\tNo. T6SS-:{absence_count}\tLog2 P/A Value:{value:.3f}\tChi-Square p-value: {p:.3f}")
