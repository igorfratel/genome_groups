#ifndef __GENOMIC_NEIGHBORHOOD_H__
#define __GENOMIC_NEIGHBORHOOD_H__

#include <iostream>
#include <string>
#include <vector>

typedef struct {
	std::string locus;
	std::string pid;
	std::string cds;
} protein_info_t;


class GenomicNeighborhood {

	private:
		std::string accession;
		std::string organism;
		std::
		std::vector<protein_info_t> seeds;
		std::vector<protein_info_t> proteins;


	public:
		typedef std::vector<protein_info_t>::iterator iterator;

		GenomicNeighborhood (std::string acession_code, std::string organism_name);

		/*Receives the locus, pid and cds of an anchor/seed protein and adds that info to the object*/
		void add_seed(std::string locus, std::string pid, std::string cds);

		/*Receives the locus, pid and cds of a protein and adds the protein to the object*/
		void add_protein(std::string locus, std::string pid, std::string cds);

		/*Returns genomic neighborhood accession code*/
		std::string get_accession();

		/*Returns name of the genomic neighborhood's organism*/
		std::string get_organism();

		/*Returns anchor/seed protein*/
		protein_info_t get_seed();

		/*Returns number of proteins in the genomic neighborhood*/
		int protein_count();

		/*Iterator for protein_info_t types, in the order that they were inserted in the object*/
		iterator begin();

		iterator end();
};

#endif
