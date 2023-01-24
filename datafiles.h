#ifndef DATAFILES_H
#define DATAFILES_H

#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include<fmt/format.h>


/*
Linear grid evaluator
*/
double linear_step(double, double, size_t, size_t);

/*
Log-spaced grid evaluator
*/
double log_step(double, double, size_t, size_t);


/*
Read the passed stream, extracting a set number of key-value pairs
into the passed map. Pairs are assumed to be listed one pair per line
with the key first and then the value, separated by whitespace.
Anything following the value will be ignored.
*/
template <typename KeyType, typename ValueType>
void read_param_entries(std::istream & in,
                        std::map<KeyType, ValueType> & params,
                        std::string & comment)
{
    std::string line;
    KeyType name;
    while (std::getline(in, line)){
        auto comment_start = line.find(comment);
        if (comment_start != std::string::npos) line.erase(comment_start);
        std::stringstream line_stream(line);
        if (line_stream >> name) line_stream >> params[name];
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

