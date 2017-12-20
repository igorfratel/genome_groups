#ifndef __PROTEIN_COLLECTION_H__
#define __PROTEIN_COLLECTION_H__

#include <vector>

//Union-Find implementation
class ProteinCollection{
    std::vector<int> parent; // parent[i] = parent of i
    std::vector<int> size;   // size[i] = number of sites in tree rooted at i
                             // Note: not necessarily correct if i is not a root node
    int count; // number of components

    public:
        void add_protein();

        void Union(int p, int q);
        int find(int p);
        bool connected(int p, int q);
        int get_count();


        ProteinCollection();

		/*Creates object with known number of proteins to be added*/
		ProteinCollection(int n);

		/*Adds protein with string identifier*/
		void add_protein(std::string node);

		/*Adds connection between two existing proteins*/
		void connect_proteins(std::string node1, std::string node2);

		/*Adds two proteins and connects them with given similarity*/
		void add_connected_proteins(std::string node1, std::string node2, double weight);

		/*Returns true if given proteins are directly connected and false otherwise*/
		bool are_connected(std::string node1, std::string node2);

		/*Returns similarity value between two existing and directly connected proteins.
		 *If not connected, returns 0.0*/
		double get_similarity(std::string node1, std::string node2);

		/*Returns number of proteins in object*/
		int get_number_proteins();

		/*Returns vector of connected components where each position is a vector of nodes in the same
		 *component*/
		std::vector<std::vector<std::string> > connected_components(double weight);

    private:


};
#endif
