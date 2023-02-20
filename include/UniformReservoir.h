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

class UniformReservoir 
{
public:
    UniformReservoir(size_t desired_sample_size); 
    ~UniformReservoir(); 

    void put(std::string element); 

    void drain(std::ofstream& sample_stream); 
    
    private:
        size_t _sample_size;
        std::string * _sample;
        int _usage;
};




