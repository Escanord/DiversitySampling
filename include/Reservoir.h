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
        Node(std::string s, double w) {
            str = s;
            weight = w;
        }

        std::string str;
        double weight;
};

struct Comp{
    bool operator()(const struct Node a, const struct Node b){
        return a.weight < b.weight;
    }
};

class Reservoir 
{
public:
    Reservoir(size_t desired_sample_size); 
    ~Reservoir(); 

    void put(std::string element, double weight); 

    void drain(std::ofstream& out); 
    
    private:
    	std::priority_queue<struct Node, std::vector<struct Node>, Comp> _pq;
        size_t _sample_size;
};




