HOW TO RUN:  

- 1: ./make  
- 2: ./neighborhood_comparer <full or partial> <args according to chosen mode>  

full --> default execution  
    -e --execution_mode full    
    -n --neighborhoods_filename  "File containing the genomic neighborhoods"   
    -s --prot_sim_filename  "File containing pairs of proteins and their similarities"  
    -f --formatted_prot_filename "File already formatted as the input for the homology detection method"  
    -p --protein_comparing  "Method for comparing proteins (default: nc)"  
    -t --prot_stringency  "Minimum similarity required to treat two proteins as a related pair (default 0)"  
    -r --neigh_stringency "Minimum threshold to display the similarity between two neighborhoods (default 0)"  
    -g --neigh_comparing  "Method for comparing genomic neighborhoods (default: porthodom method)"  
    -o --output  "Where the neighborhood similarities should be written (outputs do stdout if the filename is - or if not used)"
    -a --pairings_filename "Where the chosen pairings between proteins in the neighborhoods should be written (if not specified, does not generate a pairings file)"  


partial --> Already has the similarities between the proteins.  
    -e --execution_mode partial  
    -n --neighborhoods_filename  
    -s --prot_sim_filename  
    -l --normalize_prot_sim" "Indicates that the protein similarities file should be normalized"
    -t --prot_stringency  
    -r --neigh_stringency "Minimum threshold to display the similarity between two neighborhoods"
    -g --neigh_comparing  
    -o --output  
    -a --pairings_filename "Where the chosen pairings between proteins in the neighborhoods should be written"



Help option: -h --help  
protein scoring methods: nc or raw score in the partial mode (just use a file with the desired scores and normalize if needed).
neighborhood scoring methods: porthodom, porthodomO2, global_alignment.

prot_sim_filename format: "prot1 prot2 score" <-- whitespaces!
Output format: "accession1    cds_begin1    cds_end1    accession2    cds_begin2    cds_end2    score"  <-- tabs
Pairings format: ">accession1    cds_begin1    cds_end1    accession2    cds_begin2    cds_end2
                  prot1 prot2 sim"  <-- tabs

To add a new neighborhood scoring method: include the filename containing the scoring method functions in genome_grouping.h,
add a new "if else" clause at the genome_clustering function in the genome_grouping.cpp file comparing the genomic neighborhoods using the new scoring function. Add new files to Makefile.
