// Read a gzipped list of kmers and a list of gzipped kmer count files and define patterns for each listed kmer
// *** Assumes list and kmer files are all sorted ***
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
#include <string>

using namespace std;
using namespace myutils;

int main(const int argc, const char* argv[]) {
	if(argc!=4 && argc!=6) error("Usage: stringlist2pattern sorted-kmer-list.txt.gz kmer-sorted-gzfile-list.txt out-prefix [k min-coverage]");
    const char* kmerlistfile = argv[1];
	const char* filename = argv[2];
	const char* out_prefix = argv[3];
	// kmer size: defaults to unconstrained (0)
	const int kmerlen = (argc==4) ? 0 : atoi(argv[4]);
	if(kmerlen<0) error("The minimum kmer length is 1. k = 0 means any length");
	// Allows even kmer numbers
	// Minimum kmer count - can be
	const int mincount = (argc==4) ? 5 : atoi(argv[5]);
	if(mincount<1) error("The min-coverage must be positive");
	if(mincount>10000) error("The min-coverage cannot exceed 10000");
	// Confirm arguments (or their default values)
    if(kmerlen==0) cout << "kmer length (k) is unconstrained" << endl;
	else cout << "kmer length (k) = " << kmerlen << endl;
	cout << "Read min-coverage = " << mincount << endl;
	
	// Read the kmer file names and check the files exist
	vector<string> kmer_filename(0);
	ifstream infile(filename);
	if(!infile) {
		cout << "Could not open kmer-file-list file " << filename << endl;
		error("");
	}
	unsigned int i,j;
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
	const unsigned int n = kmer_filename.size();
    cout << "Opened input files" << endl;

    // Read the kmer list
    vector<string> kmer_list(0);
    zstr::ifstream kinfile(kmerlistfile);
    if(!kinfile) {
        cout << "Could not open gzipped kmer-list file " << kmerlistfile << endl;
        error("");
    }
    // Loop over the listed kmers
    for(i=0;!kinfile.eof();i++) {
        string next_kmer;
        kinfile >> next_kmer;
        if(kinfile.eof()) break; // This should not be necessary: eof() and operator! should not overlap - possible zstr problem?
        if(!kinfile) {
            cout << "Unexpected error reading kmer on line " << i+1 << " of " << kmerlistfile << endl;
            cout << "Read " << next_kmer << endl;
            error("");
        }
        if(kmerlen!=0) {
            if(next_kmer.length()!=kmerlen) {
                cout << "kmer on line " << i+1 << " of " << kmerlistfile << " had length " << next_kmer.length() << " when " << kmerlen << " expected" << endl;
                error("");
            }
        }
        kmer_list.push_back(next_kmer);
        //cout << "Read kmer " << next_kmer << endl;
    }
    //kinfile.close(); No explicit close, destructor invoked when deleted from the stack
    const unsigned int nkmers = kmer_list.size();
    cout << "Read " << nkmers << " kmers from " << kmerlistfile << endl;
    
    // Check the list is sorted
    if(!is_sorted(kmer_list.begin(),kmer_list.end())) {
        cout << "The input kmer list " << kmerlistfile << " is not sorted (or sorting is incompatible)" << endl;
        error("");
    }
    // Identify the smallest and largest kmer present in the kmer list
    string smallest_target = kmer_list.front();
    string largest_target = kmer_list.back();
    
	// Open the output files in readiness for writing
    // No longer necessary to write kmers, will match input list
	//zstr::ofstream kmerOut((string(out_prefix)+".kmer.txt.gz").c_str());
	// New output files for the pattern: an index per kmer and a key
	zstr::ofstream patternIndexOut((string(out_prefix)+".patternIndex.txt.gz").c_str());
	zstr::ofstream patternKeyOut((string(out_prefix)+".patternKey.txt.gz").c_str());
	
	// Create storage for the patterns
    const vector<bool> zero_pattern(n,0);
    vector< vector<bool> > kmer_pattern(nkmers,zero_pattern);
	// Previously: map<unsigned long long int, vector<bool>> kmer_pattern;

	// Open the files one at a time and construct the patterns
	for(i=0;i<n;i++) {
        // Start looking at the beginning of the target kmer list
        vector<string>::const_iterator next_target_kmer_it = kmer_list.begin();
        const vector<string>::const_iterator last_target_kmer_it = kmer_list.end();
        // Use this to ensure the input file is sorted and deduplicated
        string last_kmer;
        
		// Open the file using the zstr C++ wrapper for zlib
		zstr::ifstream gzin(kmer_filename[i].c_str());
		// Cycle through kmers, assuming no kmer appears more than once
		int nlines = 0;
		while(!gzin.eof()) {
			// Extract the kmer
            string next_kmer;
            gzin >> next_kmer;
            if(!gzin) {
                if(gzin.eof()) break; // This should not be necessary: eof() and operator! should not overlap - possible zstr problem?
                cout << "Unexpected error reading kmer on line " << nlines+1 << " of " << kmer_filename[i] << endl;
                cout << "Read " << next_kmer << endl;
                error("");
            }
            if(kmerlen!=0) {
                if(next_kmer.length()!=kmerlen) {
                    cout << "kmer on line " << nlines+1 << " of " << kmer_filename[i] << " had length " << next_kmer.length() << " when " << kmerlen << " expected" << endl;
                    error("");
                }
            }
            if(nlines>0) {
                if(next_kmer<=last_kmer) {
                    cout << "Next kmer " << next_kmer << " is smaller or equal to last kmer " << last_kmer << endl;
                    cout << "Input files must be sorted and deduplicated" << endl;
                    error("");
                }
            }
			// Extract the kmer count
            if(gzin.eof()) {
                cout << "Input file " << kmer_filename[i] << " ended unexpectedly on line " << nlines+1 << " between kmer and kmer count" << endl;
                error("");
            }
            int next_kmer_count;
            gzin >> next_kmer_count;
            if(!gzin) {
                cout << "Unexpected error reading kmer count on line " << nlines+1 << " of " << kmer_filename[i] << endl;
                cout << "Read " << next_kmer_count << endl;
                error("");
            }
            // If the kmer count exceeds the threshold, treat it as present
			if(next_kmer_count>=mincount) {
                if(next_kmer<smallest_target) {
                    // If next_kmer is smaller than smallest_target, ignore
                } else if(next_kmer>largest_target) {
                    // If next_kmer is larger than largest_target, close this file
                    break;
                } else {
                    // Otherwise, test whether the kmer is in the target list
                    const vector<string>::const_iterator next_found_kmer_it = find(next_target_kmer_it,last_target_kmer_it,next_kmer);
                    if(next_found_kmer_it==kmer_list.end()) {
                        // If there is no match, ignore
                    } else {
                        // Find the position of next_kmer in kmer_list
                        const unsigned int pos = next_found_kmer_it - kmer_list.begin();
                        // Record that the kmer was present in this individual
                        kmer_pattern[pos][i] = true;
                        // And update next_target_kmer - assumes both kmer_list and gzin are sorted
                        next_target_kmer_it = next_found_kmer_it;
                        ++next_target_kmer_it;
                    }
                }
            }
            last_kmer = next_kmer;
			++nlines;
		}
		//gzin.close(); No explicit close, destructor invoked when deleted from the stack at the close of this iteration of the loop
		cout << "Read " << nlines << " lines from file " << kmer_filename[i] << endl;
	}
    cout << "Completed reading in kmers and patterns" << endl;
	
	// Define the unique patterns and output all kmers, unique patterns key and pattern index
	unordered_map<vector<bool>,unsigned int> unique_pattern;
    vector< vector<bool> >::const_iterator kpatit;
	// Previously: map<unsigned long long int, vector<bool>>::const_iterator kpatit;
	unsigned int n_unique_patterns = 0;
	for(kpatit=kmer_pattern.begin();kpatit!=kmer_pattern.end();kpatit++) {
        // Has the current pattern (in kmer_pattern) been previously observed (in unique_pattern)?
		const unordered_map<vector<bool>,unsigned int>::const_iterator upatit = unique_pattern.find(*kpatit);
		int ipat;
		if(upatit==unique_pattern.end()) {
			// Create new pattern in memory
			ipat = n_unique_patterns;
			unique_pattern[*kpatit] = ipat;
			++n_unique_patterns;
			// Output new pattern to key file
			string scurrent_pattern(n,' ');
			for(j=0;j<n;j++) scurrent_pattern[j] = ((*kpatit)[j]==0) ? '0' : '1';
			patternKeyOut << scurrent_pattern << endl;
		} else {
			ipat = upatit->second;
		}
		// Output pattern index to file (no longer necessary to output kmers as provided as input)
		//kmerOut << hashtostring(current_kmer) << endl;
		patternIndexOut << ipat << endl;
	}
    cout << "Completed writing unique patterns and unique pattern key" << endl;
	
	// Close the output files (no explicit calls)
	//kmerOut.close();
	//patternIndexOut.close();
	//patternKeyOut.close();
	
	return 0;
}
