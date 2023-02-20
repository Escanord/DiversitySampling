#pragma once

#include <iomanip>
#include <cstdint>
#include <queue>
#include <iostream>
#include <tgmath.h> 
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>

struct Node {
    public:
        Node(std::string str, double score, double weight) {
            str = str;
            score = score;
            weight = weight;
        }

        std::string str;
        double score;
        double weight;
};

struct Comp{
    bool operator()(const struct Node a, const struct Node b){
        return a.score < b.score;
    }
};

class Reservoir 
{
public:
    Reservoir(size_t desired_sample_size); 
    ~Reservoir(); 

    void put(std::string element, double weight); 

    void drain(std::ofstream& sample_stream, std::ofstream& weight_stream); 
    
    private:
    	std::priority_queue<struct Node, std::vector<struct Node>, Comp> _pq;
        size_t _sample_size;
};




