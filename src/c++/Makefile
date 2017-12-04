CC = g++
CFLAGS = -Wall -std=c++11 -I -pg

default: projeto

projeto: projeto.cpp UndirectedEdgeWeightedGraph.o GenomicNeighborhood.o Hungarian.o\
	 			 protein_grouping.o genome_grouping.o
	$(CC) projeto.cpp $(CFLAGS) -o projeto UndirectedEdgeWeightedGraph.o GenomicNeighborhood.o\
	 															 Hungarian.o protein_grouping.o genome_grouping.o

protein_grouping: protein_grouping.cpp protein_grouping.h
	$(CC) $(CFLAGS) -c protein_grouping

genome_grouping: genome_grouping.cpp genome_grouping.h
	$(CC) $(CFLAGS) -c genome_grouping

UndirectedEdgeWeightedGraph: UndirectedEdgeWeightedGraph.cpp UndirectedEdgeWeightedGraph.h
	$(CC) $(CFLAGS) -c UndirectedEdgeWeightedGraph

GenomicNeighborhood: GenomicNeighborhood.cpp GenomicNeighborhood.h
	$(CC) $(CFLAGS) -c GenomicNeighborhood

Hungarian: Hungarian.cpp Hungarian.h
	$(CC) $(CFLAGS) -c Hungarian

clean:
	$(RM) projeto *.o *~
