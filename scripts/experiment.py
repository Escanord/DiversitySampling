import os
import subprocess
import argparse
import random
import time 
import statistics

import matplotlib.pyplot as plt
from collections import defaultdict
from Bio import SeqIO

UNCLASSIFIED_SPECIES = 0

# def run_kraken(fastq_file, seed):
#     # subprocess.call(["kraken2","--db", "/scratch1/zx22/bio/refseq214", "--threads", "12","--output",\
#     #     src_path+"out_race_"+str(size)+"_repeat_"+str(r), src_path+"race_"+str(size)+"_repeat_"+str(r)+".fastq"]) 
#     pass 

def run_diversity_sampling(fastq_file, sample_size, seed):
    output_file = f"./outputs/diverse-sample_seed={seed}_{os.path.basename(fastq_file)}"
    t0 = time.time_ns()
    subprocess.call(["./bin/diversesample", str(sample_size), "SE", fastq_file, output_file, "--seed", str(seed)]) 
    t1 = time.time_ns()
    return (output_file, output_file+".weights", t1 - t0)

def run_uniform_sampling(fastq_file, sample_size, seed):
    output_file = f"./outputs/uniform-sample_seed={seed}_{os.path.basename(fastq_file)}"
    t0 = time.time_ns()
    subprocess.call(["./bin/uniformsample", str(sample_size), "SE", fastq_file, output_file, "--seed", str(seed)]) 
    t1 = time.time_ns()
    return (output_file, t1 - t0)

#Parse commandlines
parser = argparse.ArgumentParser(
                    prog = 'RunExperiment',
                    description = 'What the program does',
                    epilog = '')

parser.add_argument('-f', '--fastq', help="fastq file pre sampling")

parser.add_argument('-a', '--sample_amount', help="post sampling fastq file", default=10000, type=int)  
parser.add_argument('-s', '--seed', help="seed to use", default=random.randint(0, 1 << 31), type=int)
parser.add_argument('-r', '--repetitions', help="number of times to repeat experiment", default=1, type=int)

parser.add_argument('-k', '--kraken', help="(Optional) kraken of original fastq file")
parser.add_argument('-v', '--verbose', action='store_true')  # on/off flag

args = parser.parse_args()

# Define verbose print function
def vprint(*x):
    if args.verbose:
        print(*x)

# Generate a seed for every repetition
vprint("Using sourceseed", args.seed)
random.seed(args.seed)
seeds = [random.randint(0, 1 << 31) for _ in range(args.repetitions)]

# Extract files
fastq_path = args.fastq
kraken_path = None

if (args.kraken == None):
    vprint("Running kraken...")
    # TODO: Execute kraken when needed.
    raise Exception("kraken execution not implemented")
    # kraken_file = run_kraken(fastq_file)
else:
    kraken_path = args.kraken
    vprint("Kraken file can be located at " + kraken_path)

# Identify each species form each sequence using kraken
id_to_species = defaultdict(lambda: UNCLASSIFIED_SPECIES)
true_proportion = defaultdict(lambda: 0)
all_diverse_estimates = defaultdict(lambda: list())
all_uniform_estimates = defaultdict(lambda: list())

# Extract species "true" proportion
numSequences = 0
with open(kraken_path) as infile:
    for line in infile:
        # Extract information from kraken line
        chunks = line.split('\t')
        classified = (chunks[0].strip() == 'C')  
        id = chunks[1]
        if classified and chunks[2] != None:
            species = int(chunks[2])
        else:
            species = UNCLASSIFIED_SPECIES
        # Map id to species
        id_to_species[id] = species
        true_proportion[species] += 1
        numSequences += 1
for species in true_proportion.keys():
    true_proportion[species] /= numSequences

for rep in range(args.repetitions):
    vprint(f"Running repitition #{rep+1}")
    seed = seeds[rep]
    vprint(f"seed={seed}")

    #Run the different sampling approaches
    vprint("Running uniform sampling...")
    (uniform_sample_path, uniform_time_elapsed) = run_uniform_sampling(fastq_path, args.sample_amount, seed)
    vprint(f" - Uniform sampling took {uniform_time_elapsed} ns")
    vprint(" - Uniform sample file can be located at " + uniform_sample_path)

    vprint("Running diversity sampling...")
    (diverse_sample_path, diverse_weights_path, diverse_time_elapsed) = run_diversity_sampling(fastq_path, args.sample_amount, seed) 
    vprint(f" - Diversity sampling took {diverse_time_elapsed} ns")
    vprint(" - Diversity sample file can be located at " + diverse_sample_path)
    vprint(" - Diversity sample Weights file can be located at " + diverse_weights_path)

    # Extract ids from samples 
    uniform_ids_list = list()
    with open(uniform_sample_path) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            uniform_ids_list.append(record.id)
    diverse_ids_list = list()
    with open(diverse_sample_path) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            diverse_ids_list.append(record.id)
    diverse_weights_list = list() 
    with open(diverse_weights_path) as infile:
        for line in infile:
            weight = float(line)
            diverse_weights_list.append(weight)

    assert(len(diverse_weights_list) == len(diverse_ids_list)) #Sanity check

    # Compute uniform estimate
    uniform_estimate = defaultdict(lambda: 0)
    total_count = 0
    for id in uniform_ids_list:
        uniform_estimate[id_to_species[id]] += 1
        total_count += 1
    for species in uniform_estimate.keys():
        uniform_estimate[species] /= total_count

    # Compute diverse estimate
    diverse_estimate = defaultdict(lambda: 0)
    total_weight = 0
    for (id, weight) in zip(diverse_ids_list, diverse_weights_list):
        diverse_estimate[id_to_species[id]] += weight
        total_weight += weight
    for species in diverse_estimate.keys():
        diverse_estimate[species] /= total_weight

    for species in true_proportion.keys():
        all_diverse_estimates[species].append(diverse_estimate[species])
        all_uniform_estimates[species].append(uniform_estimate[species])

# Organize results
rows = []
for species in true_proportion.keys():
    if (species == UNCLASSIFIED_SPECIES):
        continue #Don't plot the unclassified species
    true_pro = true_proportion[species]
    d_est = statistics.median(all_diverse_estimates[species])
    u_est = statistics.median(all_uniform_estimates[species])
    rows.append((species, true_pro, d_est, u_est))

rows.sort(key=lambda row: row[1], reverse=True)

# Print results to terminal
for row in [("Species", "Proportion", "Diverse Estimate (Median)", "Uniform Estimate (Median)")] + rows:
    print("{: >10} {: >25} {: >25} {: >25}".format(*row))

# filtered_rows = []
# Filter for only infrequent species

# LIMIT = 1e-3
# for row in rows:
#     (species, true_pro, d_est, u_est) = row
#     if true_pro < LIMIT:
#         filtered_rows.append(row)
# rows = filtered_rows

# Create plots

# ignore = 0
# rows = rows[ignore:]

limit = 0.002

x = [true_pro for (species, true_pro, d_est, u_est) in rows]
y_diverse = [d_est for (species, true_pro, d_est, u_est) in rows]
y_uniform = [u_est for (species, true_pro, d_est, u_est) in rows]

# plt.yscale("log")
# plt.xscale("log")

plt.xlim([rows[-1][1], limit])
# plt.ylim([0, limit])

plt.plot(x, y_diverse, "-b", label="Diverse Sampling",linestyle="",marker="+")
plt.plot(x, y_uniform, "-r", label="Uniform Sampling",linestyle="",marker="x")
plt.plot(x, x, color="black", label="Ideal Estimate",linestyle="dashed",marker="")

plt.legend(loc="upper left")

plt.xlabel("True Proportion")
plt.ylabel("Proportion Estimate")
plt.title('Estimates vs. Species Proportions')
plt.savefig("plot.png")

