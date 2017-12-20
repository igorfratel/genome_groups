#include "UnionFind.h"
#include <iostream>

int main() {
    UnionFind my_union_find(10);
    std::cout << "Count: " << my_union_find.get_count() << "\n";
    std::cout << "Union 0 and 1:\n";
    my_union_find.Union(0, 1);
    std::cout << "Are 0 and 1 connected?: " << my_union_find.connected(0, 1) << "\n";
    std::cout << "Are 1 and 0 connected?: " << my_union_find.connected(1, 0) << "\n";
    std::cout << "Find 0: " << my_union_find.find(0) << "\n";
    std::cout << "Find 1: " << my_union_find.find(1) << "\n";
    std::cout << "Count: " << my_union_find.get_count() << "\n";

    std::cout << "Are 0 and 9 connected?: " << my_union_find.connected(0, 9) << "\n";
    std::cout << "Are 1 and 9 connected?: " << my_union_find.connected(1, 9) << "\n";
    std::cout << "Find 9: " << my_union_find.find(9) << "\n";
    std::cout << "Union 0 and 9:\n";
    my_union_find.Union(0, 9);
    std::cout << "Are 0 and 9 connected?: " << my_union_find.connected(0, 9) << "\n";
    std::cout << "Are 1 and 9 connected?: " << my_union_find.connected(1, 9) << "\n";
    std::cout << "Find 0: " << my_union_find.find(0) << "\n";
    std::cout << "Find 9: " << my_union_find.find(9) << "\n";
    std::cout << "Count: " << my_union_find.get_count() << "\n";


}
