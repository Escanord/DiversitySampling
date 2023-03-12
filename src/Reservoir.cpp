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

void Reservoir::put(std::string element, long double weight, long double kde){
	
	long double r = (long double) rand() / RAND_MAX;
	long double score = pow(r, 1/weight);

	_pq.emplace(element, score, weight, kde);

	if (_pq.size() > _sample_size) {
		_pq.pop();
	}
}

void Reservoir::drain(std::ofstream& sample_stream, std::ofstream& weight_stream){

	while(!_pq.empty()) {
		const Node node = _pq.top();
		sample_stream << node.str;
		weight_stream << std::to_string(1/node.weight) << '\n';
		_pq.pop();
	}
}
