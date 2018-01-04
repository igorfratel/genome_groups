//Union-Find implementation

#include "ProteinCollection.h"

ProteinCollection::ProteinCollection() {}

ProteinCollection::ProteinCollection(int n) : protein_indexes(n) { //!!! resize bigger may be better
    protein_count = 0;
    count = n;
    parent.resize(n);
    size.resize(n);
    for (int i = 0; i < n; i++) {
        parent[i] = i;
        size[i] = 1;
    }
}

void ProteinCollection::add_protein(std::string p){
        if(get_number_proteins() == parent.size()){
            parent.resize(parent.size() + 1);
            size.resize(parent.size() + 1);
            parent[parent.size() - 1] = parent.size() - 1;
            size[size.size() - 1] = size.size() - 1;
        }
        protein_indexes.insert(p, get_number_proteins());
        protein_count++;
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
    int rootP = find(p);
    int rootQ = find(q);
    if (rootP == rootQ || rootP == -1 || rootQ == -1) return;

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
int ProteinCollection::find(std::string p) {
    int tmp = protein_indexes.get_value(p);
    int root = tmp;
    if (root == -1) return -1;
    while (root != parent[root])
        root = parent[root];
    while (tmp != root) {
        int new_tmp = parent[tmp];
        parent[tmp] = root;
        tmp = new_tmp;
    }
    return root;
}
bool ProteinCollection::are_connected(std::string p, std::string q){
    return find(p) == find(q);
}

int ProteinCollection::get_count(){
    return count;

}
