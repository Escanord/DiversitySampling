from Bio import SeqIO
from Bio.SeqIO import FastaIO
import random
import numpy as np
import os
import matplotlib.pyplot as plt


def random_shuffle(dir):
    records = get_records(dir)
    random.shuffle(records)
    with open(dir, "w") as output_handle:
        fasta_writer = FastaIO.FastaWriter(output_handle, wrap=None)
        fasta_writer.write_file(records)


def get_records(dir):
    records = list(SeqIO.parse(dir, "fasta"))
    return records


def get_distinct(samples):
    distinct = 0
    hash_table = {}

    for sample in samples:
        id = sample.id.split(".")[0]
        try:
            hash_table[id] = hash_table[id] + 1
        except KeyError as error:
            hash_table[id] = 1
            distinct = distinct + 1

    return distinct


def uniform_sampling(records, k):
    uniform_samples = []
    for record in records:
        flag = random.randint(1, (int)(len(records) / k))
        if flag == 1:
            uniform_samples.append(record)
    return get_distinct(uniform_samples)


def race_sampling(in_dir, out_dir, k):
    cmd = f"./bin/diversesample {k} SE {in_dir} {out_dir}"
    os.popen(cmd).read()

    records = get_records(out_dir)

    return get_distinct(records)


def plot_sample(num_samples, uniform, race):
    iters = np.arange(0, num_samples)
    plt.plot(iters, uniform, label="Uniform sampling")
    plt.plot(iters, race, label="Race Sampling")
    plt.legend(loc="upper left")
    plt.xlabel("Sampling size")
    plt.ylabel("Distinct Genes")
    plt.show(block=True)


def main():
    random_shuffle("accuracy/HiSeq_oneline.fasta")
    records = get_records("accuracy/HiSeq_oneline.fasta")
    uniform_samples = []
    race_samples = []

    for i in range(1, 11):
        uniform_ave = 0
        race_ave = 0
        for j in range(30):
            distinct = uniform_sampling(records, i)
            uniform_ave = uniform_ave + distinct
            distinct = race_sampling(
                "./accuracy/HiSeq_oneline.fasta", "./accuracy/HiSeq_out.fasta", i
            )
            race_ave = race_ave + distinct
        uniform_samples.append(round(uniform_ave / 30, 3))
        race_samples.append(round(race_ave / 30, 3))

    print(uniform_samples)
    print(race_samples)

    plot_sample(10, uniform_samples, race_samples)


if __name__ == "__main__":
    main()
