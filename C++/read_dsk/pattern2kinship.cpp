// Calculate the kinship matrix from patternKey.txt.gz and patternCounts.txt.gz
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
#include <iomanip>

using namespace std;
using namespace myutils;

int main(const int argc, const char* argv[]) {
	if(argc!=2) error("Usage: pattern2kinship prefix");
    const char* prefix = argv[1];

    // Initialize input streams and ensure files exist
    zstr::ifstream in_counts(string(prefix)+".patternCounts.txt.gz");
    zstr::ifstream in_key(string(prefix)+".patternKey.txt.gz");
    if(!in_counts) {
        cout << "Could not open patternCounts file " << prefix << ".patternCounts.txt.gz" << endl;
        error("");
    }
    if(!in_key) {
        cout << "Could not open patternKey file " << prefix << ".patternKey.txt.gz" << endl;
        error("");
    }

    // Initialize output streams to ensure they are writable
    zstr::ofstream out_kinship(string(prefix)+".kinship.txt.gz");
	
	// Concurrently read patternKey and patternCounts and construct kinship matrix
	vector< vector<double> > kinship(0);
    vector<double> zero_pattern(0);
    string s;
    unsigned int n = 0, totalCount = 0, i, j, nkey = 0;
    while(getline(in_key,s)) {
		vector<double> patternIn(zero_pattern);
		double tot = 0.0;
        if(n==0) {
            for(i=0;i<s.length();i++) {
                const unsigned char c = s[i];
                switch(c) {
                    case '0':
                        patternIn.push_back(0.0);
                        break;
                    case '1':
                        patternIn.push_back(1.0);
						tot += 1.0;
                        break;
                    default:
                        cout << "Unexpected character on line " << nkey+1 << " of " << prefix << ".patternKey.txt.gz: " << c << endl;
                        error("");
                }
            }
            n = patternIn.size();
            zero_pattern = vector<double>(n,0.0);
			// Point at which kinship matrix memory is allocated
			cout << "Allocating memory for " << n << " by " << n << " kinship matrix" << endl;
			kinship = vector< vector<double> >(n,zero_pattern);
			cout << "Allocated" << endl;
        } else {
            for(i=0;i<s.length();i++) {
                if(i==n) {
                    cout << "Expected pattern of length " << n << " at line " << nkey+1 << " of " << prefix << ".patternKey.txt.gz" << endl;
                    error("");
                }
                const unsigned char c = s[i];
                switch(c) {
                    case '0':
                        break;
                    case '1':
                        patternIn[i] = 1.0;
						tot += 1.0;
                        break;
                    default:
                        cout << "Unexpected character on line " << nkey+1 << " of " << prefix << ".patternKey.txt.gz: " << c << endl;
                        error("");
                }
            }
            if(i<n) {
                cout << "Expected pattern of length " << n << " at line " << nkey+1 << " of " << prefix << ".patternKey.txt.gz but read" << i << endl;
                error("");
            }
        }
		// Read next patternCount
		unsigned int count;
		in_counts >> count;
		if(in_counts.eof()) {
			cout << "Reached end of patternCounts file at line " << nkey+1 << " before end of patternKey file" << endl;
			error("");
		}
		if(in_counts.fail()) {
			cout << "Could not read entry from line " << nkey+1 << " of patternCounts" << endl;
			error("");
		}
		totalCount += count;
		
		// Update the kinship matrix
		for(i=0;i<n;i++) {
			patternIn[i] -= tot/(double)n;
			patternIn[i] *= sqrt((double)count);
		}
		for(i=0;i<n;i++) {
			for(j=0;j<=i;j++) {
				kinship[i][j] += patternIn[i]*patternIn[j];
			}
		}
		
		++nkey;
		if((nkey % 1000000)==0) cout << "Read " << nkey << " keys and counts totalling " << totalCount << " kmers" << endl;
    }
	cout << "Read " << nkey << " keys and counts totalling " << totalCount << " kmers" << endl;
	
	// Divide the kinship matrix by the total number of kmers
	for(i=0;i<n;i++) {
		for(j=0;j<=i;j++) {
			kinship[i][j] /= (double)totalCount;
		}
	}
	// Fill in the opposite triangle
	for(i=0;i<n;i++) {
		for(j=i+1;j<n;j++) {
			kinship[i][j] = kinship[j][i];
		}
	}
	
	// Write the kinship matrix
	out_kinship << setprecision(17);
	for(i=0;i<n;i++) {
		for(j=0;j<n;j++) {
			if(j>0) out_kinship << "\t";
			out_kinship << kinship[i][j];
		}
		out_kinship << endl;
	}
	cout << "Wrote kinship matrix" << endl;
	
    // Close the output gzip files (no explicit calls)
	
	return 0;
}
