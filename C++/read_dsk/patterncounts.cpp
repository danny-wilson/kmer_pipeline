// Use the patternIndex file to count the number of each patternKey
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
	if(argc!=2) error("Usage: patterncounts prefix");
    const char* prefix = argv[1];

    // Initialize input streams and ensure files exist
    zstr::ifstream in_index(string(prefix)+".patternIndex.txt.gz");
    if(!in_index) {
        cout << "Could not open patternIndex file " << prefix << ".patternIndex.txt.gz" << endl;
        error("");
    }

    // Initialize output streams to ensure they are writable
    zstr::ofstream out_counts(string(prefix)+".patternCounts.txt.gz");
	
    // Read original index
    // NB: indices are recorded in patternIndex.txt.gz in a zero-based system
    vector<unsigned int> counts(0);
    counts.push_back(0);
    unsigned int maxInd = 0, nInd = 0;
    while(true) {
        unsigned int ind;
        in_index >> ind;
        if(in_index.eof()) break;
        if(in_index.fail()) {
            error("Expected only integers in patternIndex.txt.gz");
        }
        if(ind>maxInd) {
            //cout << "ind = " << ind << " maxInd = " << maxInd << endl;
            // Assert that the indices be in increasing order
            if(ind!=maxInd+1) {
                cout << "Expected next new index to be " << maxInd+1 << " but read " << ind << endl;
                error("");
            }
            counts.push_back(1);
            maxInd = ind;
        } else {
            //cout << "ind = " << ind << " maxInd = " << maxInd << endl;
            ++counts[ind];
        }
        ++nInd;
        if((nInd%1000000)==0) {
            cout << "Read " << nInd << " indices across " << maxInd+1 << " keys" << endl;
        }
    }

    // Write the pattern counts
    unsigned int i;
    for(i=0;i<=maxInd;i++) {
        out_counts << counts[i] << endl;
    }
    cout << "Wrote " << nInd << " indices across " << maxInd << " keys" << endl;

    // Close the output gzip files (no explicit calls)
	
	return 0;
}
