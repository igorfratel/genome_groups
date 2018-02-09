#ifndef __GENOMIC_NEIGHBORHOOD_H__
#define __GENOMIC_NEIGHBORHOOD_H__

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>

typedef struct {
	std::string locus;
	std::string pid;
	int cds_begin;
	int cds_end;
} protein_info_t;


class GenomicNeighborhood {

	private:
		std::string accession;
		std::string organism;
		std::vector<protein_info_t> seeds;
		std::vector<protein_info_t> proteins;

		/*Receives the cds string in a format like "534..345" and splits it in the two composing numbers*/
		/*Returns a vector with two positions.*/
		std::vector<int> parse_cds(const std::string cds);


	public:
		typedef std::vector<protein_info_t>::iterator iterator;

		GenomicNeighborhood (const std::string &acession_code, const std::string &organism_name);

		/*Receives the locus, pid and cds of an anchor/seed protein and adds that info to the object*/
		void add_seed(const std::string &locus, const std::string &pid, const std::string &cds);

		/*Receives the locus, pid and cds of a protein and adds the protein to the object*/
		void add_protein(const std::string &locus, const std::string &pid, const std::string &cds);

		/*Returns genomic neighborhood accession code*/
		std::string get_accession();

		/*Returns name of the genomic neighborhood's organism*/
		std::string get_organism();

		/*Returns anchor/seed proteins*/
		std::vector<protein_info_t> get_seeds();

		/*Returns number of proteins in the genomic neighborhood*/
		int protein_count();

		/*Returns vector containing the first and last coordinates of the genomic neighborhood*/
		std::vector<int> get_cds();

		/*Returns vector containing the first and last coordinates of the genomic neighborhood
		*in a string "43434..5654" format*/
		std::string get_cds_string();

		/*Iterator for protein_info_t types, in the order that they were inserted in the object*/
		iterator begin();

		iterator end();
};

#endif
