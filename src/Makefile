CC = g++
CFLAGS = -Wall -Wextra -Wpedantic -std=c++11 -I -pg

default: neighborhood_comparator

neighborhood_comparator: neighborhood_comparator.cpp ProteinCollection.o GenomicNeighborhood.o Hungarian.o\
	 			 protein_grouping.o genome_grouping.o
	$(CC) neighborhood_comparator.cpp $(CFLAGS) -o neighborhood_comparator ProteinCollection.o GenomicNeighborhood.o\
	 															 Hungarian.o protein_grouping.o genome_grouping.o

protein_grouping: protein_grouping.cpp protein_grouping.h
	$(CC) $(CFLAGS) -c protein_grouping

genome_grouping: genome_grouping.cpp genome_grouping.h
	$(CC) $(CFLAGS) -c genome_grouping

ProteinCollection: ProteinCollection.cpp ProteinCollection.h
	$(CC) $(CFLAGS) -c ProteinCollection

GenomicNeighborhood: GenomicNeighborhood.cpp GenomicNeighborhood.h
	$(CC) $(CFLAGS) -c GenomicNeighborhood

Hungarian: Hungarian.cpp Hungarian.h
	$(CC) $(CFLAGS) -c Hungarian

clean:
	$(RM) neighborhood_comparator *.o *~
