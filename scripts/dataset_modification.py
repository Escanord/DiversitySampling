import random
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from collections import defaultdict

def get_records(dir):
    records = list(SeqIO.parse(dir, "fasta"))
    return records

def diversify_shuffle(dir):
    records = get_records(dir)
    random.shuffle(records)
    species = records[0].id.split('.')[0]
    new_records = []
    rare_species = defaultdict(lambda:[])
    
    for record in records:
        id = record.id.split('.')[0]
        if id == 'V_cholerae_HiSeq':
                new_records.append(record)
        else:
            if len(rare_species[id]) < 10:
                rare_species[id].append(record)
    
    for _, species in rare_species.items():
        for record in species:
            new_records.append(record)
    random.shuffle(new_records)

    with open("accuracy/HiSeq_rare_diversified_oneline.fasta", "w") as output_handle:
        fasta_writer = FastaIO.FastaWriter(output_handle, wrap=None)
        fasta_writer.write_file(new_records)
        
diversify_shuffle("accuracy/HiSeq_oneline.fasta")