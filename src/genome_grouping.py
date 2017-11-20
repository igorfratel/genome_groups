import subprocess
import networkx as nx
from networkx.algorithms import matching

def genome_clustering(neighborhoods_file, clusters_file, method):
    my_graph = nx.Graph()
    partition_index = 0
    if (method == "simple"):
        #PorthoDom scoring method
        with open(neighborhoods_file, 'r') as f:
            for line in f:
                line = line.split()
                for i in range (1, len(line)):
                    #sets the nodes of each partition of the graph
                    my_graph.add_node(line[i], bipartite=partition_index)
                partition_index += 1

        #Names each partition 
        left_nodes = set(n for n,d in my_graph.nodes(data=True) if d['bipartite']==0)
        right_nodes = set(n for n,d in my_graph.nodes(data=True) if d['bipartite']==1)

        #Connect each node of partition 1 to each node of partition 2
        for node1 in left_nodes:
            for node2 in right_nodes:
                key = (node1, node2)
                key_reverse = (node2, node1)
                if key in clusters_file:
                    my_graph.add_edge(node1, node2, weight=clusters_file[key])
                elif key_reverse in clusters_file:
                    my_graph.add_edge(node1, node2, weight=clusters_file[key_reverse])
                else:
                    my_graph.add_edge(node1, node2, weight=0.000001) #Default edge weight (zeroes are ignored by the matching algorithm)

        #Applies the maximum weight matching algorithm
        match = nx.matching.max_weight_matching(my_graph)

        #Applies the scoring formula
        sum_temp = 0
        for key in match:
            pair = (key, match[key])
            if pair in clusters_file:
                sum_temp += clusters_file[pair]
        sum_temp = sum_temp/2 #The matching algorithm counts each edge twice
        gen_sim = sum_temp/max(len(left_nodes), len(right_nodes)) #Contains the score between the two genomic neighborhoods

