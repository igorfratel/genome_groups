#ifndef __UTILS_H__
#define __UTILS_H__

#include <vector>
#include <string>

/**
 *Receives a string and delimiters.
 *Splits the string in a vector according to delimiters
 */
std::vector<std::string> split(const std::string& in, const std::string& delim);

#endif