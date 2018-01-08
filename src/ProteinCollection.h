#ifndef __PROTEIN_COLLECTION_H__
#define __PROTEIN_COLLECTION_H__

#include <unordered_map>

//Union-Find implementation
class ProteinCollection{
    std::unordered_map<std::string, std::string> parent; // parent[i] = parent of i
    std::unordered_map<std::string, int> size;   // size[i] = number of sites in tree rooted at i
                                                // Note: not necessarily correct if i is not a root node
    int count; // number of components

    public:
        ProteinCollection();

		/*Adds protein with string identifier*/
		void add_protein(std::string p);

		/*Adds connection between two existing proteins*/
		void connect_proteins(std::string p, std::string q);

		/*Adds two proteins and connects them.*/
		void add_connected_proteins(std::string p, std::string q);

		/*Returns true 1 given proteins are directly connected and 0 otherwise*/
		int are_connected(std::string p, std::string q);

		/*Returns similarity value between two existing and directly connected proteins.
		 *If not connected, returns 0.0*/
		double get_similarity(std::string p, std::string q);

		/*Returns number of proteins in object*/
		int get_number_proteins();

    private:
        int protein_count;
        std::string find(std::string p);
        int get_count();


};
#endif
