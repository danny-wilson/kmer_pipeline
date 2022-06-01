// Read a list of kmer count files and define patterns for each kmer
#include <stdio.h>
#include <map>
#include <iostream>
#include <stdlib.h>
#include <myerror.h>
#include <limits.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <bitset>
#include <unordered_map>
#include <array>
#include <zlib.h>
#include "zstr.hpp"

// The kmer length is a global variable shared by main() and hashtostring()
// When read in, main() mutates kmerlen and imposes a range of 1-31
unsigned int kmerlen = 31;

using namespace std;
using namespace myutils;

// Convert ACGT to 0123 to allow encoding of the 31-kmer and a long long int
char encodeACGT(const char c) {
	switch(c) {
		case 'A':
		return '0';
		case 'C':
		return '1';
		case 'G':
		return '2';
		case 'T':
		return '3';
		default:
		return '?';
	}
}

string hashtostring(const unsigned long long int x_in) {
	// The kmer is kmerlen long
	char ret[kmerlen];
	unsigned long long int x = x_in;
	const unsigned long long int four = 4;
	int i;
	for(i=kmerlen-1;i>=0;i--) {
		unsigned long long int digit = x % four;
		x /= four;
		if(digit==(unsigned long long int)0) {
			ret[i] = 'A';
		} else if(digit==(unsigned long long int)1) {
			ret[i] = 'C';
		} else if(digit==(unsigned long long int)2) {
			ret[i] = 'G';
		} else if(digit==(unsigned long long int)3) {
			ret[i] = 'T';
		} else {
			error("hashtostring(): unexpected error");
		}
	}
	return string(ret);
}

int main(const int argc, const char* argv[]) {
	if(argc!=3 && argc!=5) error("Usage: kmer2pattern kmer-file-list.txt out-prefix [k min-coverage]");
	const char* filename = argv[1];
	const char* out_prefix = argv[2];
	// Kmer size
	kmerlen = (argc==3) ? 31 : atoi(argv[3]);
	// Any longer than this and cannot use long long int for the hash table
	if(kmerlen>31) error("The maximum kmer length is 31");
	if(kmerlen<1) error("The minimum kmer length is 1");
	// Allow even kmer numbers??
	// Minimum kmer count - can be
	const int mincount = (argc==3) ? 5 : atoi(argv[4]);
	if(mincount<1) error("The min-coverage must be positive");
	if(mincount>10000) error("The min-coverage cannot exceed 10000");
	// Confirm arguments (or their default values)
	cout << "Kmer length (k) = " << kmerlen << endl;
	cout << "Read min-coverage = " << mincount << endl;
	
	// Read the kmer file names and check the files exist
	vector<string> kmer_filename(0);
	ifstream infile(filename);
	if(!infile) {
		cout << "Could not open kmer-file-list file " << filename << endl;
		error("");
	}
	int i,j;
	for(i=0;!infile.eof();i++) {
		string kmer_filename_in;
		infile >> kmer_filename_in;
		if(infile.eof()) break;
		if(!infile.good()) {
			cout << "Problem reading filename in " << filename << " line " << i+1 << endl;
			error("");
		}
		// Test whether the file exists
		ifstream kmer_infile(kmer_filename_in.c_str());
		if(!kmer_infile) {
			cout << "Could not open kmer file " << kmer_filename_in << endl;
			error("");
		}
		kmer_infile.close();
		kmer_filename.push_back(kmer_filename_in);
	}
	const int n = kmer_filename.size();
	const int MAXN = n;
	
	// Open the output files in readiness for writing
	zstr::ofstream kmerOut((string(out_prefix)+".kmer.txt.gz").c_str());
	// New output files for the pattern: an index per kmer and a key
	zstr::ofstream patternIndexOut((string(out_prefix)+".patternIndex.txt.gz").c_str());
	zstr::ofstream patternKeyOut((string(out_prefix)+".patternKey.txt.gz").c_str());
	
	// Create storage for the patterns
	map<unsigned long long int, vector<bool>> kmer_pattern;

	// Storage for reading kmers and kmer counts
	const int pbuff = 100;
	char buff[pbuff], skmer[kmerlen], skmer_count[pbuff-kmerlen];
	// Open the files one at a time and construct the patterns
	for(i=0;i<n;i++) {
		// Using the zstr C++ wrapper for zlib
		zstr::ifstream gzin(kmer_filename[i].c_str());
		// Cycle through kmers, assuming no kmer appears more than once
		int nlines = 0;
		while(gzin.getline(buff,pbuff)) {
			// Extract the kmer
			for(j=0;j<kmerlen;j++) {
				skmer[j] = encodeACGT(buff[j]);
				if(skmer[j]=='?') {
					cout << "Unexpected character '" << buff[j] << "' in file " << kmer_filename[i] << endl;
					error("");
				}
			}
			char* pend;
			const unsigned long long int next_kmer = strtoull(skmer,&pend,4);
			// Extract the kmer count
			for(;j<pbuff;j++) skmer_count[j-kmerlen] = buff[j];
			const int next_kmer_count = atoi(skmer_count);
			// Sanity check
			if(next_kmer>=LLONG_MAX) {
				cout << "The kmer " << buff << " caused an overflow" << endl;
				error("");
			}

			if(next_kmer_count>=mincount) {
				// Record that the kmer was present in this individual
				// Find the kmer if it was previously seen
				map<unsigned long long int, vector<bool>>::iterator patit = kmer_pattern.find(next_kmer);
				if(patit==kmer_pattern.end()) {
					// As the kmer was not previously been seen, create an empty pattern
					static const vector<bool> zero_pattern(MAXN,0);
					pair<unsigned long long int, vector<bool>> element(next_kmer,zero_pattern);
					patit = kmer_pattern.insert(element).first;
				}
				// Edit the presence flag for this particular individual
				patit->second[i] = 1;
			}
			
			++nlines;
		}
		//gzin.close(); No explicit close, destructor invoked when deleted from the stack at the close of this iteration of the loop
		cout << "Read " << nlines << " lines from file " << kmer_filename[i] << endl;
	}
	
	// Define the unique patterns and output all kmers, unique patterns key and pattern index
	unordered_map<vector<bool>,unsigned int> unique_pattern;
	map<unsigned long long int, vector<bool>>::const_iterator kpatit;
	unsigned int n_unique_patterns = 0;
	for(kpatit=kmer_pattern.begin();kpatit!=kmer_pattern.end();kpatit++) {
		unordered_map<vector<bool>,unsigned int>::iterator upatit;
		const unsigned long long int current_kmer = kpatit->first;
		const vector<bool> current_pattern = kpatit->second;
		upatit = unique_pattern.find(current_pattern);
		int ipat;
		if(upatit==unique_pattern.end()) {
			// Create new pattern in memory
			ipat = n_unique_patterns;
			unique_pattern[current_pattern] = ipat;
			++n_unique_patterns;
			// Output new pattern to key file
			string scurrent_pattern(n,' ');
			for(j=0;j<n;j++) scurrent_pattern[j] = (current_pattern[j]==0) ? '0' : '1';
			patternKeyOut << scurrent_pattern << endl;
		} else {
			ipat = upatit->second;
		}
		// Output kmer and pattern index to file
		kmerOut << hashtostring(current_kmer) << endl;
		patternIndexOut << ipat << endl;
	}
	
	// Close the output files (no explicit calls)
	//kmerOut.close();
	//patternIndexOut.close();
	//patternKeyOut.close();
	
	return 0;
}
