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
	ifstream in_keySizeA(string(in_prefixA)+".patternKeySize.txt");
	ifstream in_keySizeB(string(in_prefixB)+".patternKeySize.txt");
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
	if(!in_keySizeA) {
		cout << "Could not open patternKeySize file A " << in_prefixA << ".patternKeySize.txt" << endl;
		error("");
	}
	if(!in_keySizeB) {
		cout << "Could not open patternKeySize file B " << in_prefixB << ".patternKeySize.txt" << endl;
		error("");
	}

    // Initialize output streams to ensure they are writable
    zstr::ofstream out_index(string(out_prefix)+".patternIndex.txt.gz");
    zstr::ofstream out_key(string(out_prefix)+".patternKey.txt.gz");
	ofstream out_keySize(string(out_prefix)+".patternKeySize.txt");
	
	// Read key sizes: eventually nkeyA0 and nkeyB should match these sizes
	unsigned int keySizeA, keySizeB;
	in_keySizeA >> keySizeA;
	in_keySizeB >> keySizeB;
	in_keySizeA.close();
	in_keySizeB.close();
	cout << "Reading " << in_prefixA << " key of size " << keySizeA << endl;
	
    // Read all patternKey entries for A
	// Create a hash table (unordered_map) whose entries are the index in the new combined key
	// Set the number of entries to the sum of the number in the two unmerged keys: aims to avoid rehashing
	unordered_map< vector<bool>, unsigned int > keyA(keySizeA+keySizeB);
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
			vector<bool> vb_patternIn(zero_pattern);
			for(j=0;j<n;j++) if(patternIn[j]==1) vb_patternIn[j] = true;
			keyA.emplace(vb_patternIn,nkeyA);
            ++nkeyA;
        } else {
			vector<bool> vb_patternIn(zero_pattern);
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
                        vb_patternIn[j] = true;
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
			keyA.emplace(vb_patternIn,nkeyA);
            ++nkeyA;
        }
		// Write to out_key
		out_key << s << endl;
    }
	cout << "Read " << nkeyA << " keys" << endl;
	
    // Read original index A and output to new index 'C'
    while(true) {
        unsigned int ind;
        in_indexA >> ind;
        if(in_indexA.eof()) break;
        out_index << ind << endl;
    }

	cout << "Reading " << in_prefixB << " key of size " << keySizeB << endl;
    // Read all patternKey entries in B, line-by-line, and remap their indices, increasing the size of keyA on the fly
    const unsigned int nkeyA0 = nkeyA;
    unsigned int nkeyB = 0;
    vector<unsigned int> BtoA(0);
    while(getline(in_keyB,s)) {
        // First read one pattern from the key
        vector<bool> keyB(zero_pattern);
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
		const unordered_map< vector<bool>, unsigned int >::const_iterator it = keyA.find(keyB);
        unsigned int ind;
        if(it==keyA.end()) {
            // For a new pattern, add it to keyA
			keyA.emplace(keyB,nkeyA);
            // And map the B index to the (expanded) A index
            ind = nkeyA;
            ++nkeyA;
			// Write to out_key
			out_key << s << endl;
        } else {
            // For an existing pattern, map the B index to the (expanded) A index
            ind = it->second;
        }
        BtoA.push_back(ind);
    }
	cout << "Read " << nkeyB << " keys" << endl;
	cout << "Wrote " << nkeyA << " keys to " << out_prefix << endl;
	out_keySize << nkeyA << endl;
	out_keySize.close();

    // Read original index B and append to new index 'C'
    while(true) {
        unsigned int ind;
        in_indexB >> ind;
        if(in_indexB.eof()) break;
        if(ind>=BtoA.size()) error("Index B too big");
        const unsigned int newind = BtoA[ind];
        out_index << newind << endl;
    }
    
    // Close the output gzip files (no explicit calls)
	
	return 0;
}
