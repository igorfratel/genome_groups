//Union-Find implementation

#include "ProteinCollection.h"

ProteinCollection::ProteinCollection() {
    protein_count = 0;
}

void ProteinCollection::add_protein(std::string p){
    if(parent.insert(std::make_pair(p, p)).second){
        size.insert(std::make_pair(p, 1));
        protein_count++;
    }
}

void ProteinCollection::add_connected_proteins(std::string p, std::string q){
    add_protein(p);
    add_protein(q);
    connect_proteins(p, q);
}


double ProteinCollection::get_similarity(std::string p, std::string q){
    if(are_connected(p, q)) return 1.0;
    return 0.0;
}

int ProteinCollection::get_number_proteins(){
    return protein_count;
}

void ProteinCollection::connect_proteins(std::string p, std::string q) {
    std::string rootP = find(p);
    std::string rootQ = find(q);
    if (rootP == rootQ || rootP == "0" || rootQ == "0") return;

    // make smaller root point to larger one
    if (size[rootP] < size[rootQ]) {
        parent[rootP] = rootQ;
        size[rootQ] += size[rootP];
    }
    else {
        parent[rootQ] = rootP;
        size[rootP] += size[rootQ];
    }
    count--;
}
std::string ProteinCollection::find(std::string p) {
    std::string root;
    try {
        root = parent.at(p);
    }
    catch(...) {
        return "0";
    }
    while (root != parent[root])
        root = parent[root];
    while (p != root) {
        std::string new_p = parent[p];
        parent[p] = root;
        p = new_p;
    }
    return root;
}
int ProteinCollection::are_connected(std::string p, std::string q){
    return (find(p) == find(q) && find(p) != "0");
}

int ProteinCollection::get_count(){
    return count;

}
