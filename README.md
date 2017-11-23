COMO EXECUTAR:

- 1: ./make
- 2: ./python parse_neighborhood.py "arquivo com vizinhanças gênicas"
- 3: ./projeto <0 ou 1 ou 2> <args referentes ao modo escolhido>

Modo 0 --> //default execution  
           //neighborhoods_file, formatted_prot_file, protein_homology_method, num_prot, stringency, 

Modo 1 --> //Already has the similarities between the proteins. Needs to cluster them and the genomic neighborhoods.
           //neighborhoods_file, prot_sim_file, num_prot, stringency, genome_clustering_method


Modo 2 --> //Already has clustered proteins. Just needs to cluster the genomic neighborhoods
           //neighborhoods_file, prot_clusters_file, genome_clustering_method
