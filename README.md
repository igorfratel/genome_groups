Transformar o UndirectedEdgeWeightedGraph em um template pra poder usar com as vizinhanças

Usar tabelas associando strings à inteiros, ao invés de usar as strings diretamente

Passar por referência p n fazer cópia


COMO EXECUTAR:

- 1: ./make
- 2: ./python parse_neighborhood.py "arquivo com vizinhanças gênicas"
- 3: ./projeto <0 ou 1 ou 2> <args referentes ao modo escolhido>

Modo 0 --> //Already has clustered proteins. Just needs to cluster the genomic neighborhoods
           //neighborhoods_file, prot_clusters_file, genome_clustering_method

Modo 1 --> //Already has the similarities between the proteins. Needs to cluster them and the genomic neighborhoods.
           //neighborhoods_file, prot_sim_file, num_prot, stringency, genome_clustering_method

Modo 2 --> //default execution
           //neighborhoods_file, formatted_prot_file, protein_homology_method, num_prot, stringency, genome_clustering_method
