import os
import argparse

from Bio import SeqIO

#Parse commandlines
parser = argparse.ArgumentParser(
                    prog = 'RunExperiment',
                    description = 'What the program does',
                    epilog = '')

parser.add_argument('-f', '--fastq', help="fastq")
parser.add_argument('-o', '--output', help="output")
parser.add_argument('-k', '--kraken', help="(Optional) kraken of original fastq file")
parser.add_argument('-v', '--verbose', action='store_true')  # on/off flag

args = parser.parse_args()

# Define verbose print function
def vprint(*x):
    if args.verbose:
        print(*x)

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

# Extract species "true" proportion
classified_ids = set()
with open(kraken_path) as infile:
    for line in infile:
        # Extract information from kraken line
        chunks = line.split('\t')
        classified = (chunks[0].strip() == 'C')  
        id = chunks[1]
        if classified and chunks[2] != None:
            classified_ids.add(id)

classified = 0
unclassified = 0
with open(args.output, "w") as output_handle:
    with open(args.fastq, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "fastq"):
            if record.id in classified_ids:
                SeqIO.write([record], output_handle, "fastq")
                classified += 1
            else:
                unclassified += 1

print(f"Precent classified: {classified/(classified + unclassified)*100}%")

    