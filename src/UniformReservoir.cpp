#include "UniformReservoir.h"

/*
Copyright 2023, Adam Zawierucha, All rights reserved. 
Free for research use. For commercial use, contact 
Rice University Invention & Patent or the author

*/


UniformReservoir::UniformReservoir(size_t desired_sample_size){
	_sample_size = desired_sample_size;
	_usage = 0;
	_sample = (std::string *) malloc(desired_sample_size * sizeof(std::string));
}

UniformReservoir::~UniformReservoir(){
	// delete &_pq;
}

void UniformReservoir::put(std::string element){
	if (_usage < _sample_size) {
		_sample[_usage++] = element;
	} else {
		int r = rand() % _sample_size;
		_sample[r] = element;
	}
}

void UniformReservoir::drain(std::ofstream& sample_stream){
	for(int i = 0; i < _sample_size; i++) {
		sample_stream << _sample[i];
	}
}
