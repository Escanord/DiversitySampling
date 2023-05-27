# DiversitySampling


## Introduction
This repository implements a one-pass diversity sampling algorithm for DNA sequences. If you give this tool a stream of sequences or a fastq/fasta file, it will save a diverse subset of them. Sequences are saved based on their distance from sequences that are already represented in the sample. It would be very slow to explicitly calculate these distances, so we use locality-sensitive hashing (LSH) to implement the fast algorithm described in our paper. 

This is interesting if you want to quickly estimate the biodiversity of a metagenomic sample. Since different organisms often have sufficiently different sequences, our tool can ensure that each organism is represented using only a few reads. In some cases, the data reduction can be quite dramatic - in one of our experiments, RACE only needed 5% of the reads to identify 95% of the species that were present. RACE is particularly good at finding species that occur very rarely since there is a very high probability that RACE will keep an unusual or unique sequence. Most other methods cannot make this guarantee without returning very large samples. 

We tested RACE on human microbiome metagenomics data, but the tool should also work well for downsampling from other categories of labeled genetic sequences (The labels could be OTUs, species, genera, classes of bacteria, etc), provided the sequences in each category are different enough. In a nutshell, our method attempts to save a set of sequences that maximizes sequence diversity. 

## Installation
The algorithm is implemented as a command line tool written in C++. It takes in a fasta or fastq file, downsamples the file using the RACE method, and outputs a fasta or fastq file of diverse sequences. To compile the tool, clone the repo and 
```
make binaries
```
The Makefile should produce build and bin directories and output the executable file samplerace to bin/. This should work fine on most Linux systems. If something goes wrong, it is probably because your C++ compiler does not support C++11 or OpenMP. In particular, on MacOS the g++ command aliases to an outdated version of clang that does not support the -fopenmp flag. Windows does not include g++ by default, so you will need to install a compiler with OpenMP and C++11 support. 

## Algorithm and Hyperparameters
We use the RACE data structure, which is an efficient way to estimate kernel densities on streaming data. RACE is a small 2D array of integer counters indexed by a LSH function. These counters can tell whether we have already seen data that is similar to a new sequence. The key idea is that we only store sequences if we haven't seen something similar before. This gives us a diverse sample. 

### Algorithm 
RACE consists of a (B x R) array of integer counters and R different LSH functions. RACE acts as an online filter that dynamically keeps or discards sequences from a data stream. When a new sequence arrives, we hash the sequence using the R LSH functions to get R different hash values, which we use as a set of indices. We look up the integer counters at each index and take the average. The average is a *kernel density estimate*, or a measure of how many similar sequences we have seen so far. If the average is larger than a threshold, we discard the sequence since we have seen similar sequences. Otherwise, we keep the sequence. 

<img src="algorithm.png" alt="RACE Algorithm" width="800"/>

We use MinHash as the LSH function. MinHash takes two parameters - the length of the kmer to hash (k) and the number of MinHash signatures to use (n). These parameters change the sensitivity of the algorithm. If you increase (k), the Jaccard similarity between two sequences goes down. If you increase (n), RACE will become more sensitive to similarity differences. While increasing both parameters causes RACE to save more sequences, this happens for different reasons: For k, we change the interpretation of the input data while for n, we change the sensitivity of RACE to similarity differences. 

<img src="sequence_minhash.png" alt="MinHash Algorithm" width="200"/>


### Hyperparameters

There are a couple of hyperparameters that you may want to change: 

- range: This is the width (B) of the RACE array. If there are many categories or organisms that you want to sample from, increasing range might help you get more diverse results. Increasing the range is essentially free, but keeping it below 10000 may lead to faster processing times. 
- reps: This is the depth (R) of the RACE array. Increasing the reps will directly increase the time needed to process each input sequence, but you will be much less likely to accidentally discard a rare sequence. Typical values for reps are between 10 and 1000. 
- hashes: This is the number (n) of LSH functions we use for each row of the RACE array. Increasing this will directly increase the processing time but may also let you differentiate between sequences that are closer together in terms of edit distance. We recommend using only 1 hash. 
- k: This is the size (k) of each k-mer that is fed to the LSH function (MinHash). Increasing k means that we can differentiate between more similar sequences. To differentiate between species in metagenomic studies, we found that k = 16 is a good choice. If you want to differentiate between mutations or organisms within the same species, try a larger value of k. 

### Troubleshooting
If it seems like RACE isn't returning very good samples, try increasing k and increase the range. If RACE isn't returning enough samples, try increasing tau. If RACE is returning too many samples and you have already tried reducing tau, increase the reps. A more in-depth explanation of the algorithm is available in our paper. Feel free to contact the authors with any questions.

## How to run

Once you have the binaries compiled, there are two sampling methods that you can run in the /bin folder:

- bin/diversesample: executable file of RACE algorithm

- bin/uniformsample: executable file of general uniform sampling method

RACE will only require about 20 KB of RAM (in constrast to the > 10 GB needed by other diversity sampling methods such as Diginorm, coresets and buffer-based methods) and it can process about 2.5 Mbp/s on a 2016 MacBook.
```
diversesample <sample size> <format> <input> <output> [--range race_range] [--reps race_reps] [--hashes n_minhashes] [-k kmer_size]
Positional arguments: 
sample size: size of the sample to be taken from the population
format: Either PE, SE, or I for paired-end, single-end, and interleaved paired reads
input: path to input data file (.fastq or .fasta extension). For PE format, specify two files.
output: path to output sample file (same extension as input). For PE format, specify two files.
Optional arguments: 
[--range race_range]: (Optional, default 10000) Hash range for each ACE (B)
[--reps race_reps]: (Optional, default 10) Number of ACE repetitions (R)
[--hashes n_minhashes]: (Optional, default 1) Number of MinHashes for each ACE (n)
[--k kmer_size]: (Optional, default 16) Size of each MinHash k-mer (k)

Example usage:
diversesample 200 PE data/input-1.fastq data/input-2.fastq data/output-1.fastq data/output-2.fastq --range 100 --reps 50 --hashes 3 --k 5
diversesample 100 SE data/input.fastq data/output.fastq --range 100 --reps 5 --hashes 1 --k 33
diversesample 333 SE data/input.fasta data/output.fasta --range 100000 --k 20
```

We support fasta and fastq formats. For fastq files, the RACE tool supports single-end, paired-end and interleaved paired reads. RACE will decide how to parse your files based on the file extension, so be sure to input files with either the .fastq or .fasta extension. 

### Single-End Reads

You can process single-end reads using the "SE" flag. In this case, RACE needs one input file and will write one output file.

### Paired-End Reads

Process paired-end read files using the "PE" flag. In this case, RACE needs two input files and will write two output files. We require that the two input files be synchronized. That is, files 1 and 2 should contain the same number of reads in the same order. If this isn't the case, you can use [fastq-pair](https://github.com/linsalrob/fastq-pair) to produce files that RACE can read. Also, note that RACE only uses the sequences from the first paired-end read file when deciding whether to keep a read. 

### Interleaved Reads
If your reads are interleaved (i.e. the file alternately lists the forward and backward reads), then you can specify the "I" flag. RACE will decide whether to keep each pair of sequences based on the the first sequence. You only need to specify one input file and one output file. 

## Experiment & Plotting
There are two types of experiments that are implemented and you can run directly to obtain experimental plots.



### Number of Discovered Species vs. Sample Size
Experiment execution file is experiments/Diverse_Sampling_Exp.py. By executing the script, you will obtain a plot of the mean number of discovered species against the sample size.

** Note: ** To successfully execute the scripts, you have to create a Python virtual environment, or equivalent form, with required packages defined in 'requirements.txt'.

```
python experiment.py <input_file> <output_file> [max_sample_size] [repetitions]
Positional arguments: 
input_file: input file to sample from
output_file: output file to write the sample to
Optional arguments: 
[--max_sample_size]: the maximum size of the sample
[--repetitions]: number of times collecting sample from the population 

Example usage:
python experiments/Diverse_Sampling_Exp.py -i accuracy/Hiseq_oneline.fasta -o accuracy/HiSeq_out.fasta

python experiments/Diverse_Sampling_Exp.py -i accuracy/Hiseq_oneline.fasta -o accuracy/HiSeq_out.fasta -m 5 -r 5
```


### Mean Estimate Error vs. Species
Experiment execution file is script/experiment.py. By executing the script, you will obtain a plot of mean estimate error of each species' estimation against the entire population in the input file.

```
python experiment.py <fastq> [sample_amount] [seed] [repetitions]
Positional arguments: 
fastq: input file to sample from
Optional arguments: 
[--sample_amount]: the size of the sample
[--seed]: seed to be used in the sampling methods
[--repetitions]: number of times collecting sample from the population 

Example usage:
python scripts/experiment.py --fastq ./accuracy/HiSeq_rare_diversified_oneline.fasta --repetitions 50
```

## Contact 
For questions about installing or using this software, fill out an issue on GitHub and we'll do our best to help. For questions about the RACE algorithm, contact Benjamin Coleman at Rice University. If you use this software, please cite our paper. RACE is released under the MIT License. 



