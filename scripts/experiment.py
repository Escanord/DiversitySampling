import os
import subprocess
import argparse
import random
import time 
import statistics
import numpy as np
from numpy.polynomial.polynomial import polyfit

import matplotlib.pyplot as plt
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqIO import FastaIO

UNCLASSIFIED_SPECIES = 0

def stdev(data):
    if len(data) < 2:
        return 0
    return statistics.stdev(data)

# def run_kraken(fastq_file, seed):
#     # subprocess.call(["kraken2","--db", "/scratch1/zx22/bio/refseq214", "--threads", "12","--output",\
#     #     src_path+"out_race_"+str(size)+"_repeat_"+str(r), src_path+"race_"+str(size)+"_repeat_"+str(r)+".fastq"]) 
#     pass 

def run_diversity_sampling(fastq_file, sample_size, seed):
    output_file = f"./outputs/diverse-sample_seed={seed}_{os.path.basename(fastq_file)}"
    if os.path.isfile(output_file):
        print("Run diversity sampling already ran! Reusing...")
        return (output_file, output_file+".weights")
    subprocess.call(["./bin/diversesample", str(sample_size), "SE", fastq_file, output_file, "--seed", str(seed), '--range', str(29311), '--k', str(39), '--reps', str(37)]) 
    return (output_file, output_file+".weights")

def run_uniform_sampling(fastq_file, sample_size, seed):
    output_file = f"./outputs/uniform-sample_seed={seed}_{os.path.basename(fastq_file)}"
    if os.path.isfile(output_file):
        print("Run uniform sampling already ran! Reusing...")
        return (output_file)
    subprocess.call(["./bin/uniformsample", str(sample_size), "SE", fastq_file, output_file, "--seed", str(seed)]) 
    return (output_file)

#Parse commandlines
parser = argparse.ArgumentParser(
                    prog = 'RunExperiment',
                    description = 'What the program does',
                    epilog = '')

parser.add_argument('-f', '--fastq', help="fastq file pre sampling")

parser.add_argument('-a', '--sample_amount', help="post sampling fastq file", default=40, type=int)  
parser.add_argument('-s', '--seed', help="seed to use", default=random.randint(0, 1 << 31), type=int)
parser.add_argument('-r', '--repetitions', help="number of times to repeat experiment", default=1, type=int)

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

# Identify each species form each sequence using kraken
id_to_species = defaultdict(lambda: UNCLASSIFIED_SPECIES)
true_proportion = defaultdict(lambda: 0)
all_diverse_estimates = defaultdict(lambda: list())
all_uniform_estimates = defaultdict(lambda: list())
diverse_zeroes_prediction = defaultdict(lambda: 0)
uniform_zeroes_prediction = defaultdict(lambda: 0)

fastq_path = args.fastq
# Extract species "true" proportion
def get_records(dir):
    records = list(SeqIO.parse(dir, "fasta"))
    return records

def diversify_shuffle(dir):
    records = get_records(dir)
    random.shuffle(records)
    species = records[0].id.split('.')[0]
    new_records = []
    
    for record in records:
        if record.id.split('.')[0] == species:
                new_records.append(record)
        else:
            flag = random.randint(1, 10)
            if flag == 1:
                    new_records.append(record)
    random.shuffle(new_records)

    with open("accuracy/HiSeq_diversified_oneline.fasta", "w") as output_handle:
        fasta_writer = FastaIO.FastaWriter(output_handle, wrap=None)
        fasta_writer.write_file(new_records)
        
# diversify_shuffle("accuracy/HiSeq_oneline.fasta")
records = get_records("accuracy/HiSeq_rare_diversified_oneline.fasta")
numSequences = 0
for record in records:
    id = species = record.id.split('.')[0]
    # Map id to species
    id_to_species[id] = species
    true_proportion[species] += 1
    numSequences += 1
for species in true_proportion.keys():
    print(species + f' has {true_proportion[species]} sequences in the dataset')
    true_proportion[species] /= numSequences
    

for rep in range(args.repetitions):
    vprint(f"Running repitition #{rep+1}")
    seed = seeds[rep]
    vprint(f"seed={seed}")

    #Run the different sampling approaches
    vprint("Running uniform sampling...")
    (uniform_sample_path) = run_uniform_sampling(fastq_path, args.sample_amount, seed)
    # vprint(f" - Uniform sampling took {uniform_time_elapsed} ns")
    vprint(" - Uniform sample file can be located at " + uniform_sample_path)

    vprint("Running diversity sampling...")
    (diverse_sample_path, diverse_weights_path) = run_diversity_sampling(fastq_path, args.sample_amount, seed) 
    # vprint(f" - Diversity sampling took {diverse_time_elapsed} ns")
    vprint(" - Diversity sample file can be located at " + diverse_sample_path)
    vprint(" - Diversity sample Weights file can be located at " + diverse_weights_path)

    # Extract ids from samples 
    uniform_ids_list = list()
    with open(uniform_sample_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            id = record.id.split('.')[0]
            uniform_ids_list.append(id)
    diverse_ids_list = list()
    with open(diverse_sample_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            id = record.id.split('.')[0]
            diverse_ids_list.append(id)
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
        if diverse_estimate[species] == 0:
            diverse_zeroes_prediction[species] += 1
        all_uniform_estimates[species].append(uniform_estimate[species])
        if uniform_estimate[species] == 0:
            uniform_zeroes_prediction[species] += 1
    
    # print(all_diverse_estimates)
    # print(all_uniform_estimates)
    # print(true_proportion)
    
    uniform_species_detected = 0
    diverse_species_detected = 0
    for species in true_proportion.keys():
        if uniform_estimate[species] > 0:
            uniform_species_detected += 1
        if diverse_estimate[species] > 0:
            diverse_species_detected += 1 
    # print(f"Uniform species detected: {uniform_species_detected}")    
    # print(f"Diverse species detected: {diverse_species_detected}")    

# Organize results
rows = []
for species in true_proportion.keys():
    if (species == UNCLASSIFIED_SPECIES):
        continue # Don't plot the unclassified species
    true_pro = true_proportion[species]
    uniform_zeroes_prediction[species] /= args.repetitions
    diverse_zeroes_prediction[species] /= args.repetitions
    rows.append((species, true_pro, all_diverse_estimates[species], all_uniform_estimates[species]))

#Filter
# filtered_rows = []
# for row in rows:
#     (species, true_pro, d_est, u_est) = row
#     if (d_est > 0 and u_est > 0):
#         filtered_rows.append(row)
# rows = filtered_rows

# Print results to terminal
rows.sort(key=lambda row: row[1], reverse=True)

# rows = rows[100:]

# for row in [("Species", "Proportion", "Diverse Estimate (Mean)", "Uniform Estimate (Mean)")] + rows:
#     print("{: >10} {: >25} {: >25} {: >25}".format(*row))

err_uniform = defaultdict(lambda: list())
err_diverse = defaultdict(lambda: list())
est_uniform = defaultdict(lambda: list())
est_diverse = defaultdict(lambda: list())

# det_uniform = defaultdict(lambda: list())
# det_diverse = defaultdict(lambda: list())
# num_species = defaultdict(lambda: 0)
x = set()
for (species, true_pro, d, u) in rows:
    x.add(species)
    err_diverse[species] += [abs(d_est - true_pro) for d_est in d]
    err_uniform[species] += [abs(u_est - true_pro) for u_est in u]
    est_uniform[species] += u
    est_diverse[species] += d

    # num_species[true_pro] += 1
    # d_delta = 0
    # if (sum(d) > 0):
    #     d_delta = 1
    # u_delta = 0
    # if (sum(u) > 0):
    #     u_delta = 1
    # species_detected[true_pro] = (count + 1, d_detected + d_delta, u_detected + u_delta)

x = list(x)
x.sort()

# Create plots
fig, ax1 = plt.subplots()

ax2 = ax1.twinx()

plt.title('Error vs. Species Proportions')
ax1.set_xlabel("True Proportion")
ax1.set_ylabel("Mean Estimate Error (abs)")
ax2.set_ylabel("Zeroes prediction error")
plt.subplots_adjust(bottom=0.3)

zeroes_plot = ax2.errorbar([s.split('_')[1] for s in x], 
                           [proportion for _,proportion in uniform_zeroes_prediction.items()], 
                           color='tab:pink',
                           linestyle="None",
                           marker='s',
                           markeredgecolor='tab:pink',
                           markerfacecolor='None')

zeroes_plot = ax2.errorbar([s.split('_')[1] for s in x], 
                           [proportion for _,proportion in diverse_zeroes_prediction.items()], 
                           color='tab:cyan',
                           linestyle="None",
                           marker='s',
                           markeredgecolor='tab:cyan',
                           markerfacecolor='None')

#set parameters for tick labels
ax1.tick_params(axis='x', which='major', rotation = 90)

ax1.errorbar([s.split('_')[1] for s in x], 
    [statistics.mean(err_uniform[t]) for t in x],
    yerr=[stdev(err_uniform[t]) for t in x],
    capsize= 3,
    color="tab:pink", 
    label="Uniform Sampling",
    linestyle="", 
    marker="."
)

ax1.errorbar([s.split('_')[1] for s in x], 
    [statistics.mean(err_diverse[t]) for t in x],
    yerr=[stdev(err_diverse[t]) for t in x],
    capsize= 3,
    color="tab:cyan", 
    label="Diverse Sampling",
    linestyle="", 
    marker="."
)

ax1.legend(loc="upper left")
plt.savefig("error-plot.png")
plt.clf()

# Plot
plt.title('Estimate vs. Species Proportions')
plt.xlabel("True Proportion")
plt.ylabel("Mean Estimate")

plt.errorbar(x, 
    [statistics.mean(est_uniform[t]) for t in x],
    yerr=[stdev(est_uniform[t]) for t in x],
    capsize= 3,
    color="red", 
    label="Uniform Sampling",
    linestyle="", 
    marker="."
)

# b, m = polyfit(x, [statistics.mean(est_uniform[t]) for t in x], 1)
# plt.plot(x, [b + m * t for t in x ], '-', color="red", label=f"Uniform fit: {m}x + {b}")

# plt.errorbar(x, 
#     [statistics.mean(est_diverse[t]) for t in x],
#     yerr=[stdev(est_diverse[t]) for t in x],
#     capsize= 3,
#     color="blue", 
#     label="Diverse Sampling",
#     linestyle="", 
#     marker="."
# )

# b, m = polyfit(x, [statistics.mean(est_diverse[t]) for t in x], 1)
# plt.plot(x, [b + m * t for t in x ], '-', color="blue", label=f"Diverse fit: {m}x + {b}")

# plt.legend(loc="upper left")
# plt.plot(x, x, color="black", label="Ideal Estimate",linestyle="dashed",marker="")

# plt.savefig("estimate-plot.png")
# plt.clf()

# plt.title('Species Detected')
# plt.xlabel("True Proportion")
# plt.ylabel("Number of Speceis Detected")

# plt.plot(x,
#     [species_detected[t][0] for t in x],
#     color="gray",
#     label="Species Count",
#     linestyle="", 
#     marker="o"
# )

# plt.plot(x,
#     [species_detected[t][2] for t in x],
#     color="red",
#     label="Uniform Sampling",
#     linestyle="", 
#     marker="x"
# )

# plt.plot(x,
#     [species_detected[t][1] for t in x],
#     color="blue",
#     label="Diverse Sampling",
#     linestyle="", 
#     marker="+"
# )

# plt.legend(loc="upper right")
# plt.yscale("log")
# plt.savefig("detection-plot.png")
