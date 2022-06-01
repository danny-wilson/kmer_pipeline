// Read a bunch of dsk files and tabulate the kmer counts
#include <stdio.h>
#include <map>
#include <iostream>
#include <stdlib.h>
#include <myerror.h>
#include <limits.h>
#include <fstream>
#include <string>
#include <zlib.h>
#include "zstr.hpp"
#include <vector>
#include <algorithm>

using namespace std;
using namespace myutils;

int main(const int argc, const char* argv[]) {
    if(argc!=2) error("Usage: sort_strings kmers-and-counts.txt.gz");
    const char* kmerlistfile = argv[1];
    
    // Open the file using the zstr C++ wrapper for zlib
    zstr::ifstream gzin(kmerlistfile);
    if(!gzin) {
        cout << "Could not open kmers-and-counts file " << kmerlistfile << endl;
        error("");
    }
    
    vector< pair<string, int> > kmer_list(0);
    // Cycle through kmers and counts, assuming no kmer appears more than once
    int nlines = 0;
    while(!gzin.eof()) {
        // Extract the kmer
        string next_kmer;
        gzin >> next_kmer;
        if(!gzin) {
            if(gzin.eof()) break; // This should not be necessary: eof() and operator! should not overlap - possible zstr problem?
            cout << "Unexpected error reading kmer on line " << nlines+1 << " of " << kmerlistfile << endl;
            cout << "Read " << next_kmer << endl;
            error("");
        }
        // Extract the kmer count
        if(gzin.eof()) {
            cout << "Input file " << kmerlistfile << " ended unexpectedly on line " << nlines+1 << " between kmer and kmer count" << endl;
            error("");
        }
        int next_kmer_count;
        gzin >> next_kmer_count;
        if(!gzin) {
            cout << "Unexpected error reading kmer count on line " << nlines+1 << " of " << kmerlistfile << endl;
            cout << "Read " << next_kmer_count << endl;
            error("");
        }
        // Add them to the lists
        kmer_list.push_back( make_pair(next_kmer, next_kmer_count) );
        ++nlines;
    }
    //gzin.close(); No explicit close, destructor invoked when deleted from the stack at the close of this iteration of the loop
    //cout << "Read " << nlines << " lines from file " << kmerlistfile << endl;
    
    // Now sort on strings
    sort(kmer_list.begin(),kmer_list.end());
    
    // Output the total number of kmers to the standard error
    cerr << "Total\t" << nlines << endl;
    
	if(nlines>0) {
		cout << kmer_list[0].first << "\t" << kmer_list[0].second << endl;
		string *last_kmer = &(kmer_list[0].first);

		int i;
		for(i=1;i<nlines;i++) {
			// Check for duplicates
			if(kmer_list[i].first==*last_kmer) {
				cerr << "Warning: Duplicated kmer: " << last_kmer << endl;
			} else {
				cout << kmer_list[i].first << "\t" << kmer_list[i].second << endl;
				last_kmer = &(kmer_list[i].first);
			}
		}
	}
	
    return 0;
}
