import os
import subprocess
import argparse
import random
import time 

from Bio import SeqIO

UNCLASSIFIED_SPECIES = -1

# def run_kraken(fastq_file, seed):
#     # subprocess.call(["kraken2","--db", "/scratch1/zx22/bio/refseq214", "--threads", "12","--output",\
#     #     src_path+"out_race_"+str(size)+"_repeat_"+str(r), src_path+"race_"+str(size)+"_repeat_"+str(r)+".fastq"]) 
#     pass 

def run_diversity_sampling(fastq_file, seed):
    output_file = "diverse_sample_seed="+seed+"_"+fastq_file
    t0 = time.time_ns()
    subprocess.call(["./bin/diversesample", fastq_file, output_file, "--seed", seed]) 
    t1 = time.time_ns()
    return (output_file, output_file+".weights", t1 - t0)

def run_uniform_sampling(fastq_file, seed):
    output_file = "uniform_sample_seed="+seed+"_"+fastq_file
    t0 = time.time_ns()
    subprocess.call(["./bin/uniformsample", fastq_file, output_file, "--seed", seed]) 
    t1 = time.time_ns()
    return (output_file, t1 - t0)

#Parse commandlines
parser = argparse.ArgumentParser(
                    prog = 'RunExperiment',
                    description = 'What the program does',
                    epilog = '')

parser.add_argument('-f', '--fastq', help="fastq file pre sampling")

parser.add_argument('-a', '--sample_amount', help="post sampling fastq file", default=1000, type=int)  
parser.add_argument('-s', '--seed', help="seed to use", default=random.randint(0, 1 << 32), type=int)
parser.add_argument('-r', '--repetions', help="number of times to repeat experiment", default=1, type=int)

parser.add_argument('-k', '--kraken', help="(Optional) kraken of original fastq file")
parser.add_argument('-v', '--verbose', action='store_true')  # on/off flag

args = parser.parse_args()

# Define verbose print function
def vprint(*x):
    if args.verbose:
        print(*x)

# Generate a seed for every repetition
vprint("Using seed "+args.seed)
random.seed(args.seed)
seeds = [random.randint(0, 1 << 32) for _ in range(args.repetions)]

# Validate input
if (args.fastq == None and (args.sample == None and args.weights == None)):
        raise Exception("Invalid arguments. Need to provide --fastq or --sample and --weights.") 
if (args.fastq != None and (args.sample == None or args.weights == None)):
        raise Exception("Invalid arguments. Both --sample and --weights need to specified or neighter.") 
if (args.kraken == None and args.fastq == None):
    raise Exception("fastq stream must be specificied if kraken has not been specified")

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
id_to_species = dict()
species_to_true_proportion = dict()
numSequencesInKraken = 0
species_to_error = dict()

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
        if (id in ids_sampled_set):
            id_to_species[id] = species
        # Increment species count
        if (species == UNCLASSIFIED_SPECIES): #Skips unclassified species
            continue
        if (not species in species_to_true_proportion.keys()):
            species_to_true_proportion[species] = 0
        species_to_true_proportion[species] += 1
for species in species_to_true_proportion.keys():
    species_to_true_proportion[species] /= numSequencesInKraken

for rep in range(args.repetions):
    vprint(f"Running repition #{rep}")
    seed = seeds[rep]
    vprint(f"seed={seed}")

    #Run the different sampling approaches
    vprint("Running uniform sampling...")
    (uniform_sample_path, uniform_time_elapsed) = run_uniform_sampling(fastq_path, args.seed)
    vprint(f" - Uniform samploing took {uniform_time_elapsed} ns")
    vprint(" - Uniform sample file can be located at " + uniform_sample_path)

    vprint("Running diversity sampling...")
    (diverse_sample_path, diverse_weights_path, diverse_time_elapsed) = run_diversity_sampling(fastq_path, args.seed) 
    vprint(f" - Diversity sampling took {diverse_time_elapsed} ns")
    vprint(" - Diversity sample file can be located at " + diverse_sample_path)
    vprint(" - Diversity sample Weights file can be located at " + diverse_weights_path)

    ids_sampled_set = set()
    # Extract ids from samples 
    uniform_ids_list = list()
    with open(uniform_sample_path) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            uniform_ids_list.append(record.id)
            ids_sampled_set.add(record.id)
    diverse_ids_list = list()
    with open(diverse_sample_path) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            diverse_ids_list.append(record.id)
            ids_sampled_set.add(record.id)
    diverse_weights_list = list() 
    with open(diverse_weights_path) as infile:
        for line in infile:
            weight = float(line)
            diverse_weights_list.append(weight)

    assert(len(diverse_weights_list) == len(diverse_ids_list))

    # Compute uniform estimate
    species_to_uniform_estimate = dict()
    for id in uniform_ids_list:
        if not id in species_to_uniform_estimate.keys():
            species_to_uniform_estimate[id] = 0
        species_to_uniform_estimate[id] += 1
    total = len(uniform_ids_list)
    for id in uniform_ids_list:
        species_to_uniform_estimate[id] /= total
    
    # Compute diverse estimate
    species_to_diverse_estimate = dict()
    total_weight = 0
    for (id, weight) in zip(diverse_ids_list, diverse_weights_list):
        species = id_to_species[id]
        if (species == UNCLASSIFIED_SPECIES): #Skips unclassified species
            vprint("UNCLASSIFIED", id, weight)
            continue
        else:
            vprint("CLASSIFIED", id, weight)

        if (not species in species_to_diverse_estimate.keys()):
            species_to_diverse_estimate[species] = 0
        total_weight += weight
        species_to_diverse_estimate[species] += weight
    for species in species_to_diverse_estimate.keys():
        species_to_diverse_estimate[species] /= total_weight


    for species in species_to_true_proportion.keys():
        if species in species_to_error.key():
            species_to_error[species] = (0,0)
        a = species_to_true_proportion[species] or 0
        d = species_to_diverse_estimate[species] or 0
        u = species_to_uniform_estimate[species] or 0
        err = species_to_error[species]
        species_to_error[species] = (err[0] + (a - d)**2, err[1] + (a-u)**2);
        # rows.append((species, species_to_true_proportion[species] or 0, species_to_diverse_estimate[species] or 0, species_to_uniform_estimate[species] or 0))

rows = []
for species in species_to_error.keys():
    errors = species_to_error[species] or (0,0)
    rows.append((species, species_to_true_proportion[species], errors[0], errors[1]))
rows.sort(key=lambda row: row[1], reverse=True)

for row in [("Species", "Proportion", "Diverse Estimate Error", "Uniform Estimate Error")] + rows:
    print("{: >10} {: >25} {: >25} {: >25}".format(*row))