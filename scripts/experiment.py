import os
import subprocess
import argparse
from Bio import SeqIO

UNKNOWN_SPECIES = -1

def run_kraken(fastq_file):
    pass
    # subprocess.call(["kraken2","--db", "/scratch1/zx22/bio/refseq214", "--threads", "12","--output",\
    #     src_path+"out_race_"+str(size)+"_repeat_"+str(r), src_path+"race_"+str(size)+"_repeat_"+str(r)+".fastq"]) 

def run_diversity_sampling(fastq_file):
    pass



#Parse commandlines
parser = argparse.ArgumentParser(
                    prog = 'RunExperiment',
                    description = 'What the program does',
                    epilog = '')

parser.add_argument('-f', '--fastq', help="fastq file pre sampling")

parser.add_argument('-s', '--sample', help="post sampling fastq file")  
parser.add_argument('-w', '--weights', help="post sampling weights file")

parser.add_argument('-k', '--kraken', help="(Optional) kraken of original fastq file")
parser.add_argument('-v', '--verbose',
                    action='store_true')  # on/off flag

args = parser.parse_args()


def vprint(*x):
    if args.verbose:
        print(*x)

# Validate input
if (args.fastq == None and (args.sample == None and args.weights == None)):
        raise Exception("Invalid arguments. Need to provide --fastq or --sample and --weights.") 
if (args.fastq != None and (args.sample == None or args.weights == None)):
        raise Exception("Invalid arguments. Both --sample and --weights need to specified or neighter.") 
if (args.kraken == None and fastq_file == None):
    raise Exception("fastq stream must be specificied if kraken has not been specified")

# Extract files
fastq_path = args.fastq
sample_path = None
weights_path = None
kraken_path = None

if (args.sample == None):
    vprint("Diversity sampling needs to ran.")
    # TODO: Execute kraken when needed.
    raise Exception("kraken execution not implement")
    # (sample_file, weights_file) = run_diversity_sampling(fastq_file)
else:
    sample_path = args.sample 
    weights_path = args.weights
vprint("Sample file can be located at " + sample_path)
vprint("Weights file can be located at " + weights_path)

if (args.kraken == None):
    vprint("Kraken needs to be ran")
    # TODO: Execute kraken when needed.
    raise Exception("kraken execution not implemented")
    # kraken_file = run_kraken(fastq_file)
else:
    kraken_path = args.kraken
vprint("Kraken file can be located at " + kraken_path)

# Extract ids from samples 
ids_list = list()
ids_set = set()
with open(sample_path) as handle:
    for record in SeqIO.parse(handle, "fastq"):
        ids_list.append(record.id)
        ids_set.add(record.id)
vprint(ids_list)

# Extract weights
total_weight = 0
weights_list = list() 
with open(weights_path) as infile:
    for line in infile:
        weight = float(line)
        total_weight += weight
        weights_list.append(weight)

print("Number of weights:", len(weights_list))
print("Number of IDs:", len(ids_list))

assert(len(weights_list) == len(ids_list))

# Identify each speceis form each sequence using kraken
idToSpecies = dict()
speciesToProportion = dict()
numSequencesInKraken = 0

# Extract species "true" proportion
with open(kraken_path) as infile:
    for line in infile:
        numSequencesInKraken += 1
        # Extract information from kraken line
        chunks = line.split('\t')
        classified = (chunks[0] == 'C')  
        id = chunks[1]
        if (classified):
            species = int(chunks[2])
            # Map id to species
            if (id in ids_set):
                idToSpecies[id] = species
            # Increment species count
            if (not species in speciesToProportion.keys()):
                speciesToProportion[species] = 0
            speciesToProportion[species] += 1
        else: 
            if (id in ids_set):
                idToSpecies[id] = UNKNOWN_SPECIES
        # Increment species count
for species in speciesToProportion.keys():
    speciesToProportion[species] /= numSequencesInKraken

# Extract estimated species estimate
speciesToEstimate = dict()
for (id, weight) in zip(ids_list, weights_list):
    species = idToSpecies[id]
    if (not species in speciesToEstimate.keys()):
        speciesToEstimate[species] == 0
    speciesToEstimate += weight
for species in speciesToEstimate.keys():
    speciesToEstimate[species] /= total_weight

print("Species\tProportion\tEstimate")
for species in speciesToEstimate.keys():
    print(species, '\t', speciesToProportion[species], '\t', speciesToEstimate[species])
for species in speciesToProportion.keys():
    if species in speciesToEstimate.keys():
        continue
    print(species, '\t', speciesToProportion[species], '\t',  0)