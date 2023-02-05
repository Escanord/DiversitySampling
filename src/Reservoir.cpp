#include "Reservoir.h"

/*
Copyright 2023, Adam Zawierucha, All rights reserved. 
Free for research use. For commercial use, contact 
Rice University Invention & Patent or the author

*/


Reservoir::Reservoir(size_t desired_sample_size){
	_sample_size = desired_sample_size;
}

Reservoir::~Reservoir(){
	// delete &_pq;
}

void Reservoir::put(std::string element, double weight){
	
	double r = (double) rand() / RAND_MAX;
	double score = pow(r, 1/weight);

	_pq.emplace(element, score);

	if (_pq.size() > _sample_size) {
		_pq.pop();
	}
}

void Reservoir::drain(std::ofstream& out){

	while(!_pq.empty()) {
		out << _pq.top().str;
		_pq.pop();
	}
}
