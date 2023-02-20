#include "io.h"
#include "UniformReservoir.h"

#include <chrono>
#include <string>
#include <cstring>
#include <algorithm>

/*
Copyright 2023, Adam Zawierucha, All rights reserved. 
Free for research use. For commercial use, contact 
Rice University Invention & Patent or the author

*/


/*
Three types of reads: paired, interleaved and single

For single reads just do normally 
For interleaved reads just do normally but save "chunks"
For paired reads, assume they're in order (if they're not,you can use fastq-pair) 
and "rescue" saved reads


Desired interface:
samplerace tau SE input output <flags>
samplerace tau PE input1 input2 output1 output2 <flags>
samplerace tau I input output <flags>

*/

int main(int argc, char **argv){

    if (argc < 4){
        std::clog<<"Usage: "<<std::endl; 
        std::clog<<"uniformsample <sample_size> <format> <input> <output>"; 
        std::clog<<" [--seed random_seed]"<<std::endl; 
        std::clog<<"Positional arguments: "<<std::endl; 
        std::clog<<"sample_size: integer representing how many elements to sample"<<std::endl; 
        std::clog<<"format: Either PE, SE, or I for paired-end, single-end, and interleaved paired reads"<<std::endl; 
        std::clog<<"input: path to input data file (.fastq or .fasta extension). For PE format, specify two files."<<std::endl; 
        std::clog<<"output: path to output sample file (same extension as input). For PE format, specify two files."<<std::endl; 
        
        std::clog<<"Optional arguments: "<<std::endl; 
        std::clog<<"[--seed random_seed]: (Optional, default 0) The random seed to configure hash functions with"<<std::endl;

        std::clog<<std::endl<<"Example usage:"<<std::endl; 
        std::clog<<"samplerace 100 PE data/input-1.fastq data/input-2.fastq data/output-1.fastq data/output-2.fastq"<<std::endl; 
        std::clog<<"samplerace 200 SE data/input.fastq data/output.fastq -seed 100"<<std::endl; 
        std::clog<<"samplerace 300 SE data/input.fasta data/output.fasta"<<std::endl; 
        return -1; 
    }


    // POSITIONAL ARGUMENTS
    double sample_size = std::stoi(argv[1]);
    int format; // ENUM: 1 = unpaired, 2 = interleaved, 3 = paired
    if (std::strcmp("SE",argv[2]) == 0){
        format = 1;
    } else if (std::strcmp("I",argv[2]) == 0){
        format = 2; 
    } else if (std::strcmp("PE",argv[2]) == 0){
        format = 3; 
        if (argc < 7){
            std::cerr<<"For paired-end reads, please specify the input and output files as:"<<std::endl; 
            std::cerr<<"input1.fastq input2.fastq output1.fastq output2.fastq"<<std::endl; 
            return -1; 
        }
    } else {
        std::cerr<<"Invalid format, please specify either SE, PE, or I"<<std::endl; 
        return -1;
    }

    // open the correct file streams given the format
    std::ifstream datastream1;
    std::ofstream samplestream1;
    UniformReservoir reservoir1 = NULL;
    std::ifstream datastream2;
    std::ofstream samplestream2;
    UniformReservoir reservoir2 = NULL;

    if (format != 3){
        datastream1.open(argv[3]);
        samplestream1.open(argv[4]);
        reservoir1 = UniformReservoir(sample_size);
    } else {
        datastream1.open(argv[3]);
        datastream2.open(argv[4]);
        samplestream1.open(argv[5]);
        samplestream2.open(argv[6]);
        reservoir1 = UniformReservoir(sample_size);
        reservoir2 = UniformReservoir(sample_size);
    }

    // determine file extension
    std::string filename(argv[3]); 
    std::string file_extension = "";
    size_t idx = filename.rfind('.',filename.length()); 
    if (file_extension == "fq"){
        file_extension = "fastq"; 
    }
    if (idx != std::string::npos){
        file_extension = filename.substr(idx+1, filename.length() - idx); 
    } else {
        std::cerr<<"Input file does not appear to have any file extension."<<std::endl; 
        return -1; 
    }
    if (file_extension != "fasta" && file_extension != "fastq"){
        std::cerr<<"Unknown file extension: "<<file_extension<<std::endl; 
        std::cerr<<"Please specify either a file with the .fasta or .fastq extension."<<std::endl; 
        return -1; 
    }

    // OPTIONAL ARGUMENTS
    unsigned int seed = clock();

    for (int i = 0; i < argc; ++i){
         if (std::strcmp("--seed",argv[i]) == 0){
            if ((i+1) < argc){
                seed = std::stoi(argv[i+1]);
            } else {
                std::cerr<<"Invalid argument for optional parameter --seed"<<std::endl; 
                return -1;
            }
        }
    }

    srand(seed);

    // done parsing information. Begin uniform sampling algorithm algorithm: 

    // buffer for sequences and fasta/fastq chunks
    std::string sequence;
    std::string chunk1;
    std::string chunk2;

    int t = 0;

    do{
        bool success = false; 
        int c = datastream1.peek(); 
        if (c == EOF) {
            if (datastream1.eof()){
                continue; 
            }
        }

        switch(format){
            case 1: // 1 = unpaired
            success = SequenceFeaturesSE(datastream1, sequence, chunk1, file_extension);
            break; 
            case 2: // 2 = interleaved
            success = SequenceFeaturesI(datastream1, sequence, chunk1, file_extension); 
            break; 
            case 3: // 3 = paired
            success = SequenceFeaturesPE(datastream1, datastream2, sequence, chunk1, chunk2, file_extension);
            break; 
        }
        if (!success) continue;

        long double weight = 1;
  
        switch(format){
            case 1: // 1 = unpaired
            case 2: // 2 = interleaved
            reservoir1.put(chunk1);
            break; 
            case 3: // 3 = paired
            reservoir1.put(chunk1);
            reservoir2.put(chunk2);
            break; 
        }
        
    }
    while(datastream1);

    switch(format){
        case 1: // 1 = unpaired
        case 2: // 2 = interleaved
        reservoir1.drain(samplestream1);
        break; 
        case 3: // 3 = paired
        reservoir1.drain(samplestream1);
        reservoir2.drain(samplestream2);
        break; 
    }
}
