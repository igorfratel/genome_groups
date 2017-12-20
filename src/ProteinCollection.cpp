//Union-Find implementation

#include "UnionFind.h"

UnionFind::UnionFind(int n) {
    count = n;
    parent.resize(n);
    size.resize(n);
    for (int i = 0; i < n; i++) {
        parent[i] = i;
        size[i] = 1;
    }
}

void UnionFind::connect_proteins(std::string p, std::string q) {
    int rootP = find(p);
    int rootQ = find(q);
    if (rootP == rootQ) return;

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
int UnionFind::find(int p) {
    int root = p;
    while (root != parent[root])
        root = parent[root];
    while (p != root) {
        int newp = parent[p];
        parent[p] = root;
        p = newp;
    }
    return root;
}
bool UnionFind::connected(int p, int q){
    return find(p) == find(q);
}

int UnionFind::get_count(){
    return count;

}
