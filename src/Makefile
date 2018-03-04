CXX := g++
CXXFLAGS := -std=c++11 -Ofast -msse2 -march=native -Wall -Wextra -Wpedantic
LDFLAGS := -g

all: neighborhood_comparer

neighborhood_comparer: neighborhood_comparer.o ProteinCollection.o GenomicNeighborhood.o Hungarian.o\
	 			 protein_grouping.o porthodom_scoring.o porthodomO2_scoring.o global_alignment.o genome_grouping.o
	$(CXX) $^ -o $@ $(LDFLAGS)

protein_grouping.o: protein_grouping.cpp protein_grouping.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

genome_grouping.o: genome_grouping.cpp genome_grouping.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

porthodom_scoring.o: porthodom_scoring.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

porthodomO2_scoring.o: porthodomO2_scoring.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

global_alignment.o: global_alignment.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

ProteinCollection.o: ProteinCollection.cpp ProteinCollection.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

GenomicNeighborhood.o: GenomicNeighborhood.cpp GenomicNeighborhood.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

Hungarian.o: Hungarian.cpp Hungarian.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) neighborhood_comparer *.o *~
