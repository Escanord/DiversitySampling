import os
import subprocess
import argparse

repeat = 10


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
    vprint("Diversity sampling does not need to be ran.")
    sample_path = args.sample 
    weights_path = args.weights

if (args.kraken == None):
    vprint("Kraken needs to be ran")
    # TODO: Execute kraken when needed.
    raise Exception("kraken execution not implemented")
    # kraken_file = run_kraken(fastq_file)
else:
    vprint("Kraken file can be located at " + args.kraken)
    kraken_path = args.kraken

# Extract ids from samples 
ids_list = list()
with open(sample_path) as infile:
    for line in infile:
        if line.startswith('@'):
            ids_list.append(line.split(' ', 1)[0][1:])

# Extract weights
weights_list = list() 
with open(weights_path) as infile:
    for line in infile:
        weight = float(line)
        weights_list.append(weight)

assert(len(weights_list) == len(ids_list))

# Identify each speceis form each sequence using kraken
ids_set = set(ids_list)
idToSpecies = dict()
speciesToCount = dict()
numSequencesInKraken = 0
with open(kraken_path) as infile:
    for line in infile:
        # Extract information from kraken line
        chunks = line.split('\t')
        classified = (chunks[0] == 'C')
        id = chunks[1]
        species = int(chunks[2])
        # Map id to species
        idToSpecies[id] = species
        # Increment species count
        numSequencesInKraken += 1
        if (not species in speciesToCount.keys()):
            speciesToCount[species] = 0
        speciesToCount[species] += 1

