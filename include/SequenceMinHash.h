#pragma once 

#include <limits>
#include <vector>
#include <string>

#include <iostream>

#include "MurmurHash.h"


class SequenceMinHash {
private: 
	int _numhashes; 
	unsigned int _seed;
public: 
	SequenceMinHash(int number_of_hashes, unsigned int seed); 
	void getHash(size_t k, const std::string& sequence, int* hashes); 
	unsigned int internalHash(int input, int seed); 
}; 



