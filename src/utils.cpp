#include "utils.h"

/**
 *Receives a string and delimiters.
 *Splits the string in a vector according to delimiters
 */
std::vector<std::string> split(const std::string& in, const std::string& delim) {
   std::string::size_type start = in.find_first_not_of(delim), end = 0;

   std::vector<std::string> out;
   while(start != in.npos) {
      end = in.find_first_of(delim, start);
      if(end == in.npos) {
         out.push_back(in.substr(start));
         break;
      } else {
         out.push_back(in.substr(start, end-start));
      }
      start = in.find_first_not_of(delim, end);
   }
   return out;
}