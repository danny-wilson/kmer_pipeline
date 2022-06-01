// Merge patternKeys and Indicies A and B into C
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
	if(argc!=4) error("Usage: patternmerge in-prefixA in-prefixB out-prefix");
    const char* in_prefixA = argv[1];
    const char* in_prefixB = argv[2];
	const char* out_prefix = argv[3];

    // Initialize input streams and ensure files exist
    zstr::ifstream in_indexA(string(in_prefixA)+".patternIndex.txt.gz");
    zstr::ifstream in_indexB(string(in_prefixB)+".patternIndex.txt.gz");
    zstr::ifstream in_keyA(string(in_prefixA)+".patternKey.txt.gz");
    zstr::ifstream in_keyB(string(in_prefixB)+".patternKey.txt.gz");
    if(!in_indexA) {
        cout << "Could not open patternIndex file A " << in_prefixA << ".patternIndex.txt.gz" << endl;
        error("");
    }
    if(!in_indexB) {
        cout << "Could not open patternIndex file B " << in_prefixB << ".patternIndex.txt.gz" << endl;
        error("");
    }
    if(!in_keyA) {
        cout << "Could not open patternKey file A " << in_prefixA << ".patternKey.txt.gz" << endl;
        error("");
    }
    if(!in_keyB) {
        cout << "Could not open patternKey file B " << in_prefixB << ".patternKey.txt.gz" << endl;
        error("");
    }

    // Initialize output streams to ensure they are writable
    zstr::ofstream out_index(string(out_prefix)+".patternIndex.txt.gz");
    zstr::ofstream out_key(string(out_prefix)+".patternKey.txt.gz");
    
    // Read all patternKey entries for A
    vector< vector<bool> > keyA(0);
    vector<bool> zero_pattern(0);
    string s;
    unsigned int n = 0, nkeyA = 0, i, j;
    while(getline(in_keyA,s)) {
        if(n==0) {
            vector<unsigned short int> patternIn(0);
            for(i=0;i<s.length();i++) {
                const unsigned char c = s[i];
                switch(c) {
                    case '0':
                        patternIn.push_back(0);
                        break;
                    case '1':
                        patternIn.push_back(1);
                        break;
                    default:
                        cout << "Unexpected character on line " << nkeyA+1 << " of " << in_prefixA << ".patternKey.txt.gz: " << c << endl;
                        error("");
                }
            }
            n = patternIn.size();
            zero_pattern = vector<bool>(n,false);
            keyA.push_back(zero_pattern);
            for(j=0;j<n;j++) if(patternIn[j]==1) keyA[nkeyA][j] = true;
            ++nkeyA;
        } else {
            keyA.push_back(zero_pattern);
            for(i=0,j=0;i<s.length();i++) {
                if(j==n) {
                    cout << "Expected pattern of length " << n << " at line " << nkeyA+1 << " of " << in_prefixA << ".patternKey.txt.gz" << endl;
                    error("");
                }
                const unsigned char c = s[i];
                switch(c) {
                    case '0':
                        ++j;
                        break;
                    case '1':
                        keyA[nkeyA][j] = true;
                        ++j;
                        break;
                    default:
                        cout << "Unexpected character on line " << nkeyA+1 << " of " << in_prefixA << ".patternKey.txt.gz: " << c << endl;
                        error("");
                }
            }
            if(j<n) {
                cout << "Expected pattern of length " << n << " at line " << nkeyA+1 << " of " << in_prefixA << ".patternKey.txt.gz but read" << j << endl;
                error("");
            }
            ++nkeyA;
        }
    }
    
    // Read original index A and output to new index 'C'
    while(true) {
        unsigned int ind;
        in_indexA >> ind;
        if(in_indexA.eof()) break;
        out_index << ind << endl;
    }

    // Read all patternKey entries in B, line-by-line, and remap their indices, increasing the size of keyA on the fly
    const unsigned int nkeyA0 = nkeyA;
    unsigned int nkeyB = 0;
    vector<unsigned int> BtoA(0);
    while(getline(in_keyB,s)) {
        // First read one pattern from the key
        vector<bool> keyB = zero_pattern;
        for(i=0,j=0;i<s.length();i++) {
            if(j==n) {
                cout << "Expected pattern of length " << n << " at line " << nkeyB+1 << " of " << in_prefixB << ".patternKey.txt.gz" << endl;
                error("");
            }
            const unsigned char c = s[i];
            switch(c) {
                case '0':
                    ++j;
                    break;
                case '1':
                    keyB[j] = true;
                    ++j;
                    break;
                default:
                    cout << "Unexpected character on line " << nkeyB+1 << " of " << in_prefixB << ".patternKey.txt.gz: " << c << endl;
                    error("");
            }
        }
        if(j<n) {
            cout << "Expected pattern of length " << n << " at line " << nkeyB+1 << " of " << in_prefixB << ".patternKey.txt.gz but read" << j << endl;
            error("");
        }
        ++nkeyB;
        
        // Look it up in keyA
        const vector< vector<bool> >::const_iterator it = find(keyA.begin(),keyA.end(),keyB);
        unsigned int ind;
        if(it==keyA.end()) {
            // For a new pattern, add it to keyA
            keyA.push_back(keyB);
            // And map the B index to the (expanded) A index
            ind = nkeyA;
            ++nkeyA;
        } else {
            // For an existing pattern, map the B index to the (expanded) A index
            ind = it - keyA.begin();
        }
        BtoA.push_back(ind);
    }
    
    // Output expanded key A to new key 'C'
    for(i=0;i<nkeyA;i++) {
        for(j=0;j<n;j++) {
            const unsigned char c = (keyA[i][j]) ? '1' : '0';
            out_key << c;
        }
        out_key << endl;
    }

    // Read original index B and append to new index 'C'
    while(true) {
        unsigned int ind;
        in_indexB >> ind;
        if(in_indexB.eof()) break;
        if(ind>=BtoA.size()) error("Index B too big");
        const unsigned int newind = BtoA[ind];
        out_index << newind << endl;
    }
    
    // Close the output files (no explicit calls)
	
	return 0;
}
