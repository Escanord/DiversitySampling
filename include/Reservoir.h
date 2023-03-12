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
        Node(std::string s, long double c, long double w, long double k) {
            str = s;
            score = c;
            weight = w;
            kde = k;
        }

        std::string str;
        long double score;
        long double weight;
        long double kde;
};

struct Comp{
    public:
        bool operator()(const struct Node a, const struct Node b){
            return a.score > b.score;
        }
};

class Reservoir 
{
public:
    Reservoir(size_t desired_sample_size); 
    ~Reservoir(); 

    void put(std::string element, long double weight, long double kde); 

    void drain(std::ofstream& sample_stream, std::ofstream& weight_stream); 
    
    private:
    	std::priority_queue<struct Node, std::vector<struct Node>, Comp> _pq;
        size_t _sample_size;
};




