INCLUDE=-I. -Izstr -I../myutils
CPPFLAGS=-O3 -g $(INCLUDE) -std=c++0x
LDFLAGS=-O3 -g

all : kmerlist2pattern pattern2kinship patterncounts patternmerge sort_strings stringlist2count stringlist2pattern

unused : gwas_kmer kmer_info pd gwas_kmer_pattern gwas_kmer_pattern_variable_mincount gwas_kmer_pattern_15 sort_dsk sort_dsk_15 kmer2pattern

gwas_kmer : gwas_kmer.o
	g++ -o gwas_kmer gwas_kmer.o $(LDFLAGS)

kmer_info : kmer_info.o
	g++ -o kmer_info kmer_info.o $(LDFLAGS)
	
pd : pd.o
	g++ -o pd pd.o $(LDFLAGS)

gwas_kmer_pattern : gwas_kmer_pattern.o
	g++ -o gwas_kmer_pattern gwas_kmer_pattern.o $(LDFLAGS)

gwas_kmer_pattern_variable_mincount : gwas_kmer_pattern_variable_mincount.o
	g++ -o gwas_kmer_pattern_variable_mincount gwas_kmer_pattern_variable_mincount.o $(LDFLAGS)

gwas_kmer_pattern_15 : gwas_kmer_pattern_15.o
	g++ -o gwas_kmer_pattern_15 gwas_kmer_pattern_15.o $(LDFLAGS)

sort_dsk_15 : sort_dsk_15.o
	g++ -o sort_dsk_15 sort_dsk_15.o $(LDFLAGS)

sort_dsk : sort_dsk.o
	g++ -o sort_dsk sort_dsk.o $(LDFLAGS)

kmer2pattern : kmer2pattern.o
	g++ -o kmer2pattern kmer2pattern.o $(LDFLAGS) -lz

kmerlist2pattern : kmerlist2pattern.o
	g++ -o kmerlist2pattern kmerlist2pattern.o $(LDFLAGS) -lz

patternmerge : patternmerge.o
	g++ -o patternmerge patternmerge.o $(LDFLAGS) -lz

patterncounts : patterncounts.o
	g++ -o patterncounts patterncounts.o $(LDFLAGS) -lz

pattern2kinship : pattern2kinship.o
	g++ -o pattern2kinship pattern2kinship.o $(LDFLAGS) -lz

sort_strings : sort_strings.o
	g++ -o sort_strings sort_strings.o $(LDFLAGS) -lz

stringlist2pattern : stringlist2pattern.o
	g++ -o stringlist2pattern stringlist2pattern.o $(LDFLAGS) -lz

stringlist2count : stringlist2count.o
	g++ -o stringlist2count stringlist2count.o $(LDFLAGS) -lz



gwas_kmer.o : main.cpp
	g++ -c -o gwas_kmer.o main.cpp $(CPPFLAGS)

kmer_info.o : kmer_info.cpp
	g++ -c -o kmer_info.o kmer_info.cpp $(CPPFLAGS)

pd.o : pd.cpp
	g++ -c -o pd.o pd.cpp $(CPPFLAGS)

gwas_kmer_pattern.o : pattern.cpp
	g++ -c -o gwas_kmer_pattern.o pattern.cpp $(CPPFLAGS)

gwas_kmer_pattern_variable_mincount.o : pattern_variable_mincount.cpp
	g++ -c -o gwas_kmer_pattern_variable_mincount.o pattern_variable_mincount.cpp $(CPPFLAGS)

gwas_kmer_pattern_15.o : pattern15.cpp
	g++ -c -o gwas_kmer_pattern_15.o pattern15.cpp $(CPPFLAGS)

sort_dsk_15.o : sort_dsk15.cpp
	g++ -c -o sort_dsk_15.o sort_dsk15.cpp $(CPPFLAGS)

sort_dsk.o : sort_dsk.cpp
	g++ -c -o sort_dsk.o sort_dsk.cpp $(CPPFLAGS)

kmer2pattern.o : kmer2pattern.cpp
	g++ -c -o kmer2pattern.o kmer2pattern.cpp $(CPPFLAGS)

kmerlist2pattern.o : kmerlist2pattern.cpp
	g++ -c -o kmerlist2pattern.o kmerlist2pattern.cpp $(CPPFLAGS)

patternmerge.o : patternmerge.cpp
	g++ -c -o patternmerge.o patternmerge.cpp $(CPPFLAGS)

patterncounts.o : patterncounts.cpp
	g++ -c -o patterncounts.o patterncounts.cpp $(CPPFLAGS)

pattern2kinship.o : pattern2kinship.cpp
	g++ -c -o pattern2kinship.o pattern2kinship.cpp $(CPPFLAGS)

sort_strings.o : sort_strings.cpp
	g++ -c -o sort_strings.o sort_strings.cpp $(CPPFLAGS)

stringlist2pattern.o : stringlist2pattern.cpp
	g++ -c -o stringlist2pattern.o stringlist2pattern.cpp $(CPPFLAGS)

stringlist2count.o : stringlist2count.cpp
	g++ -c -o stringlist2count.o stringlist2count.cpp $(CPPFLAGS)

.PHONY : cleanobj

cleanobj :
	-rm *.o
