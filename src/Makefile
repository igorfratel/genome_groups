CC = g++
CFLAGS = -Wall -Wextra -Wpedantic -std=c++11 -I -pg

default: neighborhood_comparer

neighborhood_comparer: neighborhood_comparer.cpp ProteinCollection.o GenomicNeighborhood.o Hungarian.o\
	 			 protein_grouping.o porthodom_scoring.o porthodomO2_scoring.o genome_grouping.o
	$(CC) neighborhood_comparer.cpp $(CFLAGS) -o neighborhood_comparer ProteinCollection.o GenomicNeighborhood.o\
	 															 Hungarian.o protein_grouping.o porthodom_scoring.o\
																 porthodomO2_scoring.o genome_grouping.o

protein_grouping: protein_grouping.cpp protein_grouping.h
	$(CC) $(CFLAGS) -c protein_grouping

genome_grouping: genome_grouping.cpp genome_grouping.h
	$(CC) $(CFLAGS) -c genome_grouping

porthodom_scoring: porthodom_scoring.cpp
	$(CC) $(CFLAGS) -c porthodom_scoring

porthodomO2_scoring: porthodomO2_scoring.cpp
	$(CC) $(CFLAGS) -c porthodomO2_scoring

ProteinCollection: ProteinCollection.cpp ProteinCollection.h
	$(CC) $(CFLAGS) -c ProteinCollection

GenomicNeighborhood: GenomicNeighborhood.cpp GenomicNeighborhood.h
	$(CC) $(CFLAGS) -c GenomicNeighborhood

Hungarian: Hungarian.cpp Hungarian.h
	$(CC) $(CFLAGS) -c Hungarian

clean:
	$(RM) neighborhood_comparer *.o *~
