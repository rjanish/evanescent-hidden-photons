#ifndef DATAFILES_H
#define DATAFILES_H

#include <string>
#include <iostream>
#include <fstream>
#include <map>

#include<fmt/format.h>


/*
Linear grid evaluator
*/
double linear_step(double start, double end, size_t N, size_t index)

/*
Log-spaced grid evaluator
*/
double log_step(double start, double end, size_t N, size_t index)


/*
Read the passed stream, extracting a set number of key-value pairs
into the passed map. Pairs are assumed to be listed sequentially
and whitespace separated.
*/
template <typename KeyType, typename ValueType>
void read_param_entries(std::istream & in,
                        std::map<KeyType, ValueType> & params,
                        size_t num_entries)
{
    KeyType name;
    ValueType value;
    for (size_t counter = 0; counter < num_entries; ++counter){
        in >> name >> value;
        params[name] = value;
    }
}


/*
Print key-value pairs into a passed output stream.
Pairs are printed whitespace separated and one pair to a line.
*/
template <typename KeyType, typename ValueType>
std::ostream & write_map(std::ostream & out,
                         std::map<KeyType, ValueType> & map)
{
    typename std::map<KeyType, ValueType>::iterator it = map.begin();
    while (it != map.end()){
        out << fmt::format("{}    {}\n", it -> first, it -> second);
        ++it;
    }
    return out;
}

#endif

