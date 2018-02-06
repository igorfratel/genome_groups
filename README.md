HOW TO RUN:  

- 1: ./make  
- 2: ./python parse_neighborhood.py "file with neighborhoods"  
- 3: ./neighborhood_comparer <full or partial> <args according to chosen mode>  


ARGUMENTS FOR EACH EXECUTION MODE:  

full --> default execution  
    -e --execution_mode full  
    -n --neighborhoods_filename  
    -s --prot_sim_filename  
    -f --formatted_prot_filename  
    -p --protein_comparing  
    -m --num_prot  
    -t --stringency  
    -g --genome_comparing  
    -o --output  

partial --> Already has the similarities between the proteins.  
    -e --execution_mode partial  
    -n --neighborhoods_filename  
    -s --prot_sim_filename  
    -m --num_prot  
    -t --stringency  
    -g --genome_comparing  
    -o --output  

protein scoring methods: nc.
neighborhood scoring methods: porthodom, porthodomO2.

To add new a neighborhood scoring method: include the file containing the scoring function in genome_grouping.h,
add a new "if else" clause at the genome_clustering function in the genome_grouping.cpp file comparing the genomic neighborhoods using the new scoring function. Add new files to Makefile.
