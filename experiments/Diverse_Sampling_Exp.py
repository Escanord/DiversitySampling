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
        if flag == 1 and (len(uniform_samples) < k):
            uniform_samples.append(record)
    return get_distinct(uniform_samples)


def race_sampling(in_dir, out_dir, k):
    cmd = f"./bin/diversesample {k} SE {in_dir} {out_dir} -range 15137 --k 39 --reps 23"
    os.popen(cmd).read()

    records = get_records(out_dir)

    return get_distinct(records)


def plot_sample(num_samples, uniform, race, uniform_errors, race_errors):
    iters = np.arange(0, num_samples)
    plt.xlabel("Sampling size")
    plt.ylabel("Distinct Genes")
    plt.plot(iters, uniform, label="Uniform sampling", color = 'black')
    plt.fill_between(iters, uniform_errors[0], uniform_errors[1], edgecolor='#CC4F1B', facecolor='#FF9848')
    # plt.errorbar(iters, uniform, yerr = uniform_errors, label="Uniform sampling", color = 'blue')
    # plt.show(block=True)
    plt.plot(iters, race, label="Race Sampling", color = 'white')
    plt.fill_between(iters, race_errors[0], race_errors[1], edgecolor='#1B2ACC', facecolor='#089FFF')
    # plt.errorbar(iters, race, yerr = race_errors, label="Race Sampling", color = 'green')
    plt.legend(loc="upper left")
    plt.show(block=True)

def box_plot_sample(num_samples, uniform, race):
    # uniform plot
    iters = np.arange(0, num_samples)
    plt.xlabel("Sampling size")
    plt.ylabel("Distinct Genes")
    plt.boxplot(uniform)
    plt.errorbar(iters + 1, np.mean(uniform, axis=1), yerr=np.std(uniform, axis=1), label="Uniform sampling")
    plt.title('Uniform Sampling')
    plt.show(block=True)
    
    #race plot
    iters = np.arange(0, num_samples)
    plt.xlabel("Sampling size")
    plt.ylabel("Distinct Genes")
    plt.boxplot(race)
    plt.errorbar(iters + 1, np.mean(race, axis=1), yerr=np.std(race, axis=1), label="Race sampling")
    plt.title('Race Sampling')
    plt.show(block=True)

def main():
    random_shuffle("accuracy/HiSeq_oneline.fasta")
    records = get_records("accuracy/HiSeq_oneline.fasta")
    uniform_samples = []
    race_samples = []
    uniform_errors = [[],[]]
    race_errors = [[],[]]
    box_uniform = []
    box_race = []

    for i in range(1, 51):
        uniform_ave = 0
        uniform_min = race_min = i + 1
        uniform_max = race_max = 0
        race_ave = 0
        box_uniform.append([])
        box_race.append([])
        for j in range(20):
            distinct = uniform_sampling(records, i)
            uniform_ave = uniform_ave + distinct
            uniform_min = min(uniform_min, distinct)
            uniform_max = max(uniform_max, distinct)
            box_uniform[i - 1].append(distinct)
            
            distinct = race_sampling(
                "./accuracy/HiSeq_oneline.fasta", "./accuracy/HiSeq_out.fasta", i
            )
            race_ave = race_ave + distinct
            race_min = min(race_min, distinct)
            race_max = max(race_max, distinct)
            box_race[i - 1].append(distinct)
        
        uniform_samples.append(round(uniform_ave / 20, 3))
        race_samples.append(round(race_ave / 20, 3))
        uniform_errors[0].append(uniform_min)
        uniform_errors[1].append(uniform_max)
        race_errors[0].append(race_min)
        race_errors[1].append(race_max)

    print(uniform_samples)
    print(race_samples)

    plot_sample(50, uniform_samples, race_samples, uniform_errors, race_errors)
    # box_plot_sample(30, box_uniform, box_race)


if __name__ == "__main__":
    # a = [3.4, 6.1, 7.9, 11.0, 12.5, 13.3, 14.4, 17.1, 17.5, 19.2, 18.1, 19.8, 21.9, 21.8, 22.9, 24.3, 25.3, 25.1, 24.3, 25.5, 25.4, 27.5, 27.2, 27.3, 27.6, 27.6, 27.9, 28.1, 27.7, 28.5]
    # b = [3.0, 5.6, 8.2, 10.6, 12.5, 13.3, 16.7, 17.4, 18.5, 19.4, 20.5, 21.7, 21.7, 22.9, 24.8, 24.5, 25.6, 26.4, 25.6, 26.7, 26.8, 27.4, 27.6, 27.4, 28.1, 27.5, 28.2, 28.5, 29.1, 28.3]
    # print ([round(x/3, 3) for x in a])
    # print([round(x/3, 3) for x in b])
    main()
