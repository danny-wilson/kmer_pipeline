// Read a list of kmer count files and case control status and test for association
#include <stdio.h>
#include <map>
#include <iostream>
#include <stdlib.h>
#include <myerror.h>
#include <limits.h>
#include <fstream>
#include <vector>
#include <cmath>

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
  if(argc!=3) error("Usage: gwas_kmer kmer-pheno.txt out-prefix");
  const char* filename = argv[1];
  const char* out_prefix = argv[2];
  // Kmer size
  const int kmerlen = 31; // Any longer than this and cannot use long long int for the hash table
  // Minimum kmer count
  const int mincount = 5;

  // Read all the kmer files and the case-control status
  vector<string> kmer_filename(0);
  vector<int> kmer_phenotype(0);
  int ncontrol=0, ncase=0;
  ifstream infile(filename);
  if(!infile) {
    cout << "Could not open kmer-pheno file " << filename << endl;
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
  if((ncase+ncontrol)==0) error("No kmers or phenotypes were read");
  if(ncase==0) error("There were no cases");
  if(ncontrol==0) error("There were no controls");
  const int n = kmer_filename.size();
  if(kmer_phenotype.size()!=n) error("Difference in length of kmer files and phenotypes");

  // Open the pipes in readiness for reading
  vector<FILE*> pin(n,NULL);
  for(i=0;i<n;i++) {
    string command = "zcat ";
    command += string(kmer_filename[i]);
    pin[i] = popen(command.c_str(),"r");
    if(pin[i]==NULL) {
      cout << "Could not open stream for file " << kmer_filename[i] << endl;
      error("");
    }
  }

  // Open the output files in readiness for writing
  ofstream kmerOut((string(out_prefix)+".kmer.txt").c_str());
  ofstream nPresentCtrlOut((string(out_prefix)+".nPresentCtrl.txt").c_str());
  ofstream nPresentCaseOut((string(out_prefix)+".nPresentCase.txt").c_str());
  ofstream meanDepthCtrlOut((string(out_prefix)+".meanDepthCtrl.txt").c_str());
  ofstream meanDepthCaseOut((string(out_prefix)+".meanDepthCase.txt").c_str());
  ofstream chisqStatOut((string(out_prefix)+".chisqStat.txt").c_str());
  ofstream trendStatOut((string(out_prefix)+".trendStat.txt").c_str());
  ofstream rgrssStatOut((string(out_prefix)+".rgrssStat.txt").c_str());
  ofstream chisqLogPOut((string(out_prefix)+".chisqLogP.txt").c_str());
  ofstream trendLogPOut((string(out_prefix)+".trendLogP.txt").c_str());
  ofstream rgrssLogPOut((string(out_prefix)+".rgrssLogP.txt").c_str());

  // Cycle through files, starting with the smallest kmer (these are sorted in the input files)
  const int pbuff = 100;
  char buff[pbuff];
  char skmer[kmerlen], skmer_count[pbuff-kmerlen];
  vector<unsigned long long int> next_kmer(n);
  vector<int> next_kmer_count(n);
  unsigned long long int smallest_kmer;
  char* FGETS = NULL;
  for(i=0;i<n;i++) {
    // Read line
    if(feof(pin[i])) {
      cout << "Unexpectedly reached end of file " << kmer_filename[i] << endl;
      error("");
    }
    FGETS = fgets(buff,pbuff,pin[i]);
    if(FGETS==NULL) {
      cout << "Unexpectedly reached end of file " << kmer_filename[i] << endl;
      error("");
    }
    for(j=0;j<kmerlen;j++) {
      skmer[j] = encodeACGT(buff[j]);
	  if(skmer[j]=='?') {
          cout << "Unexpected character '" << buff[j] << "' in file " << kmer_filename[i] << endl;
          error("");
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
      cout << "The kmer " << buff << " caused an overflow" << endl;
      error("");
    }
  }

  // Cycle through all remaining kmers
  unsigned long long int nkmers = 1;
  // Kmer count per genome
  vector<int> current_kmer_count(n);
  for(;;nkmers++) {
    for(i=0;i<n;i++) {
      if(next_kmer[i]==smallest_kmer && next_kmer_count[i]>=mincount) {
	// The kmer is present in the observed count
	current_kmer_count[i] = next_kmer_count[i];
      } else {
	// Assume that if the kmer is not recorded or fall below the minimum count it is absent
	current_kmer_count[i] = 0;
      }
    }
    // Calculate the test statistics
    int nPresentCtrl = 0, nPresentCase = 0;
    double depthCtrl = 0.0, depthCase = 0.0;
    for(i=0;i<n;i++) {
      if(current_kmer_count[i]>0) {
	if(kmer_phenotype[i]==0) {
	  ++nPresentCtrl;
	  depthCtrl += (double)current_kmer_count[i];
	} else {
	  ++nPresentCase;
	  depthCase += (double)current_kmer_count[i];
	}
      }
    }
    depthCtrl /= (double)ncontrol;
    depthCase /= (double)ncase;
    // Don't bother with the tests or output if the kmer is absent or fixed in the population
    // NB:- NEED TO REVISE THIS CONDITION WHEN EXTENDING TO THE OTHER NON-PEARSON TESTS
    if(!(nPresentCtrl==0 && nPresentCase==0) && !(nPresentCtrl==ncontrol && nPresentCase==ncase)) {
      // Pearson chi-squared test
      const int nAbsentCtrl = ncontrol-nPresentCtrl;
      const int nAbsentCase = ncase-nPresentCase;
      const double expdPresentCtrl = (double)(nPresentCtrl+nPresentCase)*(double)ncontrol/(double)n;
      const double expdPresentCase = (double)(nPresentCtrl+nPresentCase)*(double)ncase/(double)n;
      const double expdAbsentCtrl = (double)(n-nPresentCtrl-nPresentCase)*(double)ncontrol/(double)n;
      const double expdAbsentCase = (double)(n-nPresentCtrl-nPresentCase)*(double)ncase/(double)n;
      const double chisqStat = pow((double)nPresentCtrl-expdPresentCtrl,2.0)/expdPresentCtrl +
	pow((double)nPresentCase-expdPresentCase,2.0)/expdPresentCase +
	pow((double)nAbsentCtrl-expdAbsentCtrl,2.0)/expdAbsentCtrl +
	pow((double)nAbsentCase-expdAbsentCase,2.0)/expdAbsentCase;
      // Calculate the p-values
      // Output to files
      kmerOut << hashtostring(smallest_kmer) << endl;
      nPresentCtrlOut << nPresentCtrl << endl;
      nPresentCaseOut << nPresentCase << endl;
      meanDepthCtrlOut << depthCtrl << endl;
      meanDepthCaseOut << depthCase << endl;
      chisqStatOut << chisqStat << endl;
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
            cout << "Unexpected character '" << buff[j] << "' in file " << kmer_filename[i] << endl;
            error("");
          }
	    }
	    for(;j<pbuff;j++) skmer_count[j-kmerlen] = buff[j];
	    char* pend;
	    // Record next kmer and kmer count
	    next_kmer[i] = strtoull(skmer,&pend,4);
	    next_kmer_count[i] = atoi(skmer_count);
	    // Sanity check
	    if(next_kmer[i]>=LLONG_MAX) {
	      cout << "The kmer " << buff << " caused an overflow" << endl;
	      error("");
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
  }

  // Close the output files
  kmerOut.close();
  nPresentCtrlOut.close();
  nPresentCaseOut.close();
  meanDepthCtrlOut.close();
  meanDepthCaseOut.close();
  chisqStatOut.close();
  trendStatOut.close();
  rgrssStatOut.close();
  chisqLogPOut.close();
  trendLogPOut.close();
  rgrssLogPOut.close();

  // Close the pipes
  for(i=0;i<n;i++) pclose(pin[i]);

  // Regurgitate to standard output the ordered map
  //  map<unsigned long long int,double>::iterator it;
  //for(it=hash.begin();it!=hash.end();it++) {
  // cout << hashtostring(it->first) << "\t" << it->second << endl;
  //}

  return 0;
}
