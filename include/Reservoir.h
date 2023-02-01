#pragma once

#include <iostream>
#include <iomanip>
#include <cstdint>
#include <queue>

using namespace std;

struct Node {
    public:
        Node(void * p, double w) {
            ptr = p;
            weight = w;
        }

        void * ptr;
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

    void process(void * element, double weight); 
    void reset(); 

    void pprint(std::ostream& out, int width = 3, bool format = true); 
    
    private:
    	priority_queue< struct Node, vector<struct Node>, Comp> _pq;
        size_t _sample_size;
};




