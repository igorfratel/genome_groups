#ifndef __PROTEIN_COLLECTION_H__
#define __PROTEIN_COLLECTION_H__

#include <vector>
#include "HashTable.h"

//Union-Find implementation
class ProteinCollection{
    HashTable protein_indexes;
    std::vector<int> parent; // parent[i] = parent of i
    std::vector<int> size;   // size[i] = number of sites in tree rooted at i
                             // Note: not necessarily correct if i is not a root node
    int count; // number of components

    public:
        ProteinCollection();

		/*Creates object with known number of proteins to be added*/
		ProteinCollection(int n);

		/*Adds protein with string identifier*/
		void add_protein(std::string p);

		/*Adds connection between two existing proteins*/
		void connect_proteins(std::string p, std::string q);

		/*Adds two proteins and connects them.*/
		void add_connected_proteins(std::string p, std::string q);

		/*Returns true if given proteins are directly connected and false otherwise*/
		bool are_connected(std::string p, std::string q);

		/*Returns similarity value between two existing and directly connected proteins.
		 *If not connected, returns 0.0*/
		double get_similarity(std::string p, std::string q);

		/*Returns number of proteins in object*/
		int get_number_proteins();

		/*Returns vector of connected components where each position is a vector of nodes in the same
		 *component*/
		std::vector<std::vector<std::string> > connected_components(double weight);

    private:
        int protein_count;
        int find(std::string p);
        int get_count();


};
#endif
