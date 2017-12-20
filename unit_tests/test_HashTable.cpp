#include "HashTable.h"
#include <iostream>

int main() {
    HashTable my_hash_table(10);
    std::string s1 = "EAA26069.1";
    std::string s2 = "EAA26070.1";
    std::string s3 = "EAA26071.1";
    std::string s4 = "EAA26072.1";
    std::string s5 = "EAL53902.1";
    std::string s6 = "EAL53903.1";
    std::string s7 = "EAL53904.1";
    std::string s8 = "EAL53905.1";
    std::string s9 = "EAL53906.1";
    my_hash_table.insert(s1, 1);
    my_hash_table.insert(s2, 2);
    my_hash_table.insert(s3, 3);
    my_hash_table.insert(s4, 4);
    my_hash_table.insert(s5, 5);
    my_hash_table.insert(s6, 6);
    my_hash_table.insert(s7, 7);
    my_hash_table.insert(s8, 8);
    my_hash_table.insert(s9, 9);
    std::cout << s1 << " " << my_hash_table.get_hash(s1) << " " << my_hash_table.get_value(s1) << "\n";
    std::cout << s2 << " " << my_hash_table.get_hash(s2) << " " << my_hash_table.get_value(s2) << "\n";
    std::cout << s3 << " " << my_hash_table.get_hash(s3) << " " << my_hash_table.get_value(s3) << "\n";
    std::cout << s4 << " " << my_hash_table.get_hash(s4) << " " << my_hash_table.get_value(s4) << "\n";
    std::cout << s5 << " " << my_hash_table.get_hash(s5) << " " << my_hash_table.get_value(s5) << "\n";
    std::cout << s6 << " " << my_hash_table.get_hash(s6) << " " << my_hash_table.get_value(s6) << "\n";
    std::cout << s7 << " " << my_hash_table.get_hash(s7) << " " << my_hash_table.get_value(s7) << "\n";
    std::cout << s8 << " " << my_hash_table.get_hash(s8) << " " << my_hash_table.get_value(s8) << "\n";
    std::cout << s9 << " " << my_hash_table.get_hash(s9) << " " << my_hash_table.get_value(s9) << "\n";
}
