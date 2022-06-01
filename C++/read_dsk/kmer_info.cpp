// Read a list of kmer count files for certain kmers, output information
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
	// The kmer is 31 long
	char ret[31];
	unsigned long long int x = x_in;
	const unsigned long long int four = 4;
	int i;
	for(i=30;i>=0;i--) {
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
	if(argc!=4) error("Usage: kmer_info kmer-pheno.txt kmer-list.txt out-prefix");
	const char* kmer_phenotype_filename = argv[1];
	const char* kmer_list_filename = argv[2];
	const char* out_prefix = argv[3];
	// Kmer size
	const int kmerlen = 31; // Any longer than this and cannot use long long int for the hash table
	// Minimum kmer count
	const int mincount = 5;
	
	// Read all the kmer files and the case-control status
	vector<string> kmer_filename(0);
	vector<int> kmer_phenotype(0);
	int ncontrol=0, ncase=0;
	ifstream infile(kmer_phenotype_filename);
	if(!infile) {
		cout << "Could not open kmer-pheno file " << kmer_phenotype_filename << endl;
	}
	int i,j;
	for(i=0;!infile.eof();i++) {
		string kmer_filename_in, kmer_phenotype_in;
		infile >> kmer_filename_in;
		if(infile.eof()) break;
		infile >> kmer_phenotype_in;
		// Test whether the file exists
		ifstream kmer_infile(kmer_filename_in.c_str());
		if(!kmer_infile) {
			cout << "Could not open kmer file " << kmer_filename_in << endl;
			error("");
		}
		kmer_infile.close();
		kmer_filename.push_back(kmer_filename_in);
		// Test whether the case control status is 0 or 1
		const int kmer_phenotype_int = atoi(kmer_phenotype_in.c_str());
		if(kmer_phenotype_int!=0 && kmer_phenotype_int!=1) {
			cout << "Phenotype (" << kmer_phenotype_int << ") was not encoded as control (0) or case (1)" << endl;
			error("");
		}
		if(kmer_phenotype_int==0) {
			++ncontrol;
		} else {
			++ncase;
		}
		kmer_phenotype.push_back(kmer_phenotype_int);
	}
	infile.close();
	if((ncase+ncontrol)==0) error("No kmers or phenotypes were read");
	const int n = kmer_filename.size();
	if(kmer_phenotype.size()!=n) error("Difference in length of kmer files and phenotypes");
	
	// Read all the kmers to process
	vector<unsigned long long int> kmer_list(0);
	const int pbuff = 100;
	char buff[pbuff], kmer[kmerlen], kmer_temp[kmerlen];
	// Could reinstate the search for the smallest and largest kmer in the list for more efficient, 
	// sorted processing of the list later rather than using std::find() to test the presence of each
	// kmer in the list, which will be inefficient if the list is very long
	//  unsigned long long int smallest_kmer, largest_kmer;
	FILE* fin = fopen(kmer_list_filename,"r");
	char* FGETS = NULL;
	if(fin==NULL) {
		cout << "Error: Could not open stream" << kmer_list_filename << endl;
		exit(13);
	}
	for(i=0;!feof(fin);i++) {
		FGETS = fgets(buff,pbuff,fin);
		if(FGETS==NULL) break;
		for(j=0;j<kmerlen;j++) {
			kmer_temp[j] = buff[j];
			kmer[j] = encodeACGT(buff[j]);
			if(kmer[j]=='?') {
				cout << "Error: Unexpected character '" << buff[j] << "' in file " << kmer_list_filename << endl;
				exit(13);
			}
		}
		// Record in hash table
		char* pend;
		unsigned long long int skmer = strtoull(kmer,&pend,4);
		//    if(i==0) {
		//      smallest_kmer = largest_kmer = skmer;
		//    } else {
		//      if(skmer<smallest_kmer) {
		//        smallest_kmer = skmer;
		//      } else if(skmer>largest_kmer) {
		//        largest_kmer = skmer;
		//      }
		//    }
		if(skmer>=LLONG_MAX) {
			cout << "The kmer " << kmer << " caused an overflow" << endl;
			error("");
		}
		kmer_list.push_back(skmer);
	}
	fclose(fin);
	if(kmer_list.size()==0) error("No kmers were read from kmer-list file");
	// Sort the kmer list
	sort(kmer_list.begin(),kmer_list.end());
	// Next smallest kmer in list
	vector<unsigned long long int>::iterator next_smallest_target_kmer = kmer_list.begin();
	
	// Open the pipes in readiness for reading
	vector<FILE*> pin(n,NULL);
	for(i=0;i<n;i++) {
		string command = "zcat ";
		command += string(kmer_filename[i]);
		pin[i] = popen(command.c_str(),"r");
		if(pin[i]==NULL) {
			cout << "Error: Could not open stream for file " << kmer_filename[i] << endl;
			exit(13);
		}
	}
	
	// Open the output files in readiness for writing
	ofstream kmer_bip((string(out_prefix)+".kmer_count.txt").c_str());
	// Output headers to bip file
	kmer_bip << "kmer";
	for(i=0;i<n;i++) {
		kmer_bip << '\t' << kmer_filename[i];
	}
	kmer_bip << endl;
	
	// Cycle through files, starting with the smallest kmer (these are sorted in the input files)
	char skmer[kmerlen], skmer_count[pbuff-kmerlen];
	vector<unsigned long long int> next_kmer(n);
	vector<int> next_kmer_count(n);
	unsigned long long int smallest_kmer;
	for(i=0;i<n;i++) {
		// Read line
		if(feof(pin[i])) {
			cout << "Error: Unexpectedly reached end of file " << kmer_filename[i] << endl;
			exit(13);
		}
		FGETS = fgets(buff,pbuff,pin[i]);
		if(FGETS==NULL) {
			cout << "Error: Unexpectedly reached end of file " << kmer_filename[i] << endl;
			exit(13);
		}
		for(j=0;j<kmerlen;j++) {
			skmer[j] = encodeACGT(buff[j]);
			if(skmer[j]=='?') {
				cout << "Error: Unexpected character '" << buff[j] << "' in file " << kmer_filename[i] << endl;
				exit(13);
			}
		}
		for(;j<pbuff;j++) skmer_count[j-kmerlen] = buff[j];
		char* pend;
		// Record next kmer and kmer count
		next_kmer[i] = strtoull(skmer,&pend,4);
		next_kmer_count[i] = atoi(skmer_count);
		// Determine smallest kmer
		if(i==0) {
			smallest_kmer = next_kmer[i];
		} else if(next_kmer[i]<smallest_kmer) {
			smallest_kmer = next_kmer[i];
		}
		// Sanity check
		if(next_kmer[i]>=LLONG_MAX) {
			cout << "Error: The kmer " << buff << " caused an overflow" << endl;
			exit(13);
		}
	}
	
	// Cycle through all remaining kmers
	unsigned long long int nkmers = 1;
	// Kmer count per genome
	vector<int> current_kmer_count(n);
	for(;;nkmers++) {
		// Calculate the test statistics and output: only if the kmer is in the list
		if(smallest_kmer==*next_smallest_target_kmer) {
			for(i=0;i<n;i++) {
				if(next_kmer[i]==smallest_kmer && next_kmer_count[i]>=mincount) {
					// The kmer is present in the observed count
					current_kmer_count[i] = next_kmer_count[i];
				} else {
					// Assume that if the kmer is not recorded or fall below the minimum count it is absent
					current_kmer_count[i] = 0;
				}
				
				if(i==0) kmer_bip << hashtostring(smallest_kmer);
				kmer_bip << '\t' << current_kmer_count[i];
				if(i==(n-1)) kmer_bip << endl;
			}
			// Advance the next target kmer
			++next_smallest_target_kmer;
		}
		
		// Iterate to next kmer
		for(i=0;i<n;i++) {
			// Advance all those files for which the next kmer was the current one
			if(next_kmer[i]==smallest_kmer) {
				// Test for end of file
				if(feof(pin[i])) {
					next_kmer[i] = 0ULL;
					next_kmer_count[i] = 0;
				} else {
					FGETS = fgets(buff,pbuff,pin[i]);
					if(FGETS==NULL) {
						next_kmer[i] = 0ULL;
						next_kmer_count[i] = 0;
					} else {
						for(j=0;j<kmerlen;j++) {
							skmer[j] = encodeACGT(buff[j]);
							if(skmer[j]=='?') {
								cout << "Error: Unexpected character '" << buff[j] << "' in file " << kmer_filename[i] << endl;
								exit(13);
							}
						}
						for(;j<pbuff;j++) skmer_count[j-kmerlen] = buff[j];
						char* pend;
						// Record next kmer and kmer count
						next_kmer[i] = strtoull(skmer,&pend,4);
						next_kmer_count[i] = atoi(skmer_count);
						// Sanity check
						if(next_kmer[i]>=LLONG_MAX) {
							cout << "Error: The kmer " << buff << " caused an overflow" << endl;
							exit(13);
						}
					}
				}
			}
		}
				
		// Update next kmer
		smallest_kmer = 0ULL; // Safe to use 0 from now on as indicator of NA
		for(i=0;i<n;i++) {
			// Determine the next smallest kmer
			if(smallest_kmer==0ULL) {
				smallest_kmer = next_kmer[i];
			} else if(next_kmer[i]!=0ULL && next_kmer[i]<smallest_kmer) {
				smallest_kmer = next_kmer[i];
			}
		}
		// If the smallest kmer had the special value 0ULL, then all files have been exhausted
		if(smallest_kmer==0ULL) break;
		
		// If the target list is exhausted, break
		if(next_smallest_target_kmer==kmer_list.end()) break;
	}
	
	// Close the output files
	kmer_bip.close();
	
	// Close the pipes
	for(i=0;i<n;i++) pclose(pin[i]);
	
	return 0;
}
