#include "HashTable.h"

HashTable::HashTable(int n) {
    table.resize(n);
    for(int i = 0; i < n; i++){
        table[i].resize(1);
        table[i][0].second = -1;
    }
}

void HashTable::insert(std::string key, int value) {
    int index = hash_function(key)%table.size();
    unsigned int j = 0;
    for(; j < table[index].size(); j++){
        if(table[index][j].first == key) return;
        if(table[index][j].second == -1){
            table[index][j].first = key;
            table[index][j].second = value;
            return;
        }
    }
    table[index].resize(table.size() + 1);
    table[index][j].first = key;
    table[index][j].second = value;
}

int HashTable::get_value(std::string key) {
    int index = hash_function(key)%table.size();
    for(unsigned int j = 0; j < table[index].size(); j++){
        if(table[index][j].first == key)
            return table[index][j].second;
    }
    return -1;
}

int HashTable::get_hash(std::string key){
    return hash_function(key)%table.size();
}
