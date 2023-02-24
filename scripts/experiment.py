import os
import subprocess
import argparse
import random
import time 

from PyGnuplot import gp
from collections import defaultdict

from Bio import SeqIO

UNCLASSIFIED_SPECIES = -1

# def run_kraken(fastq_file, seed):
#     # subprocess.call(["kraken2","--db", "/scratch1/zx22/bio/refseq214", "--threads", "12","--output",\
#     #     src_path+"out_race_"+str(size)+"_repeat_"+str(r), src_path+"race_"+str(size)+"_repeat_"+str(r)+".fastq"]) 
#     pass 

def run_diversity_sampling(fastq_file, sample_size, seed):
    output_file = f"diverse-sample_seed={seed}_{os.path.basename(fastq_file)}"
    t0 = time.time_ns()
    subprocess.call(["./bin/diversesample", str(sample_size), "SE", fastq_file, output_file, "--seed", str(seed)]) 
    t1 = time.time_ns()
    return (output_file, output_file+".weights", t1 - t0)

def run_uniform_sampling(fastq_file, sample_size, seed):
    output_file = f"uniform-sample_seed={seed}_{os.path.basename(fastq_file)}"
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
parser.add_argument('-r', '--repetions', help="number of times to repeat experiment", default=1, type=int)

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
seeds = [random.randint(0, 1 << 31) for _ in range(args.repetions)]

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
species_to_true_proportion = defaultdict(lambda: 0)
species_to_error = defaultdict(lambda: (0,0))

numSequencesInKraken = 0

# Extract species "true" proportion
with open(kraken_path) as infile:
    for line in infile:
        numSequencesInKraken += 1
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
        # Increment species count
        if (species == UNCLASSIFIED_SPECIES): #Skips unclassified species
            continue
        species_to_true_proportion[species] += 1
for species in species_to_true_proportion.keys():
    species_to_true_proportion[species] /= numSequencesInKraken

for rep in range(args.repetions):
    vprint(f"Running repitition #{rep+1}")
    seed = seeds[rep]
    vprint(f"seed={seed}")

    #Run the different sampling approaches
    vprint("Running uniform sampling...")
    (uniform_sample_path, uniform_time_elapsed) = run_uniform_sampling(fastq_path, args.sample_amount, args.seed)
    vprint(f" - Uniform sampling took {uniform_time_elapsed} ns")
    vprint(" - Uniform sample file can be located at " + uniform_sample_path)

    vprint("Running diversity sampling...")
    (diverse_sample_path, diverse_weights_path, diverse_time_elapsed) = run_diversity_sampling(fastq_path, args.sample_amount,args.seed) 
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

    assert(len(diverse_weights_list) == len(diverse_ids_list))

    # Compute uniform estimate
    species_to_uniform_estimate = defaultdict(lambda: 0)
    for id in uniform_ids_list:
        species_to_uniform_estimate[id] += 1
    total = len(uniform_ids_list)
    for id in uniform_ids_list:
        species_to_uniform_estimate[id] /= total
    
    # Compute diverse estimate
    species_to_diverse_estimate = defaultdict(lambda: 0)
    total_weight = 0
    for (id, weight) in zip(diverse_ids_list, diverse_weights_list):
        species = id_to_species[id]
        if (species == UNCLASSIFIED_SPECIES): #Skips unclassified species
            # vprint("UNCLASSIFIED", id, weight)
            continue
        # else:
        #     # vprint("CLASSIFIED", id, weight)
        total_weight += weight
        species_to_diverse_estimate[species] += weight
    for species in species_to_diverse_estimate.keys():
        species_to_diverse_estimate[species] /= total_weight


    for species in species_to_true_proportion.keys():
        a = species_to_true_proportion[species]
        d = species_to_diverse_estimate[species]
        u = species_to_uniform_estimate[species]
        err = species_to_error[species]
        species_to_error[species] = (err[0] + abs(a - d), err[1] + abs((a-u)));

# Organize results
rows = []
for species in species_to_error.keys():
    errors = species_to_error[species]
    rows.append((species, species_to_true_proportion[species], errors[0]/args.sample_amount, errors[1]/args.sample_amount))
rows.sort(key=lambda row: row[1], reverse=True)

# Print results to terminal
for row in [("Species", "Proportion", "Diverse Estimate Error", "Uniform Estimate Error")] + rows:
    print("{: >10} {: >25} {: >25} {: >25}".format(*row))

# Create plots
f1 = gp()
x = [i + 1 for i in range(len(rows))]
f1.plot([x, [r[2] for r in rows] ], "Diverse Estimate Error")
f1.plot([x, [r[3] for r in rows] ], "Uniform Estimate Error")
f1.p('myfigure.png')

