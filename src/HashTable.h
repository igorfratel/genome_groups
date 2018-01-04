#ifndef __HASH_TABLE_H__
#define __HASH_TABLE_H__

#include <string>
#include <vector>

class HashTable {
    std::vector<std::vector<std::pair<std::string,int> > > table;
    std::hash<std::string> hash_function;
    public:
        HashTable();
        HashTable(int n);
        void insert(std::string key, int value);
        int get_value(std::string key);
        int get_hash(std::string key);
};
#endif
