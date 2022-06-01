// Read a bunch of dsk files and tabulate the kmer counts
#include <stdio.h>
#include <map>
#include <iostream>
#include <stdlib.h>
#include <myerror.h>
#include <limits.h>
#include <fstream>

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
    cout << "Unexpected character: " << c << endl;
    error("");
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
  string command = "/dipro/mmm/gorm/v3/ana/Unkn/dsk/dsk-1.6066/parse_results";
  if(argc==3) {
    command = string(argv[2]);
  } else if(argc!=2) {
    error("Usage: sort_dsk filename.solid_kmers_binary [command]");
  }
  const char* filename = argv[1];
  // Kmer size
  const int kmerlen = 31; // Any longer than this and cannot use long long int for the hash table
  // Minimum kmer count
  const int mincount = 2;

  // Check the file exists
  ifstream infile(filename);
  if(!infile) {
    cout << "Could not open file " << filename << endl;
    error("");
  }
  infile.close();
  command += string(" ");
  command += string(filename);
  //FILE* pin = popen("/dipro/mmm/gorm/v3/ana/Unkn/dsk/dsk-1.6066/parse_results /dipro/mmm/gorm/v2/ana/Saur/AbxGWAS/dsk/C00001166_R00000022.solid_kmers_binary","r");
  FILE* pin = popen(command.c_str(),"r");
  if(pin==NULL) {
    cout << "Could not open stream" << endl;
    return 13;
  }

  // Read lines
  map<unsigned long long int,double> hash;
  const int pbuff = 100;
  char buff[pbuff];

  int i,j;
  char kmer[kmerlen], kmer_temp[kmerlen];
  char skmercount[pbuff-kmerlen];
  unsigned long long int smallest_kmer, largest_kmer;
  unsigned long long int nkmers = 0;
  for(i=0;!feof(pin);i++) {
    fgets(buff,pbuff,pin);
    for(j=0;j<kmerlen;j++) {
      kmer_temp[j] = buff[j];
      kmer[j] = encodeACGT(buff[j]);
    }
    for(;j<pbuff;j++) skmercount[j-kmerlen] = buff[j];
    const int kmercount = atoi(skmercount);
    // Record in hash table
    char* pend;
    unsigned long long int skmer = strtoull(kmer,&pend,4);
    if(i==0) {
      smallest_kmer = largest_kmer = skmer;
    } else {
      if(skmer<smallest_kmer) {
	smallest_kmer = skmer;
      } else if(skmer>largest_kmer) {
	largest_kmer = skmer;
      }
    }
    if(skmer>=LLONG_MAX) {
      cout << "The kmer " << kmer << " caused an overflow" << endl;
      error("");
    }
    hash[skmer] = (double)kmercount;
    nkmers += kmercount;
    if(false) {
      if((i % 1000000)==0) {
	cout << "Line " << i << " records kmer " << kmer_temp << " or " << skmer << " " << kmercount << " times (" << hashtostring(skmer) << ")" << endl;
      }
    }
  }
  pclose(pin);

  if(false) {
    cout << "Smallest kmer was " << hashtostring(smallest_kmer) << "(" << smallest_kmer << ")" << endl;
    cout << "Largest kmer was  " << hashtostring(largest_kmer) << "(" << largest_kmer << ")" << endl;
  }

  // Output the total number of kmers to the standard error
  cerr << "Total\t" << nkmers << endl;
  // Regurgitate to standard output the ordered map
  map<unsigned long long int,double>::iterator it;
  for(it=hash.begin();it!=hash.end();it++) {
   cout << hashtostring(it->first) << "\t" << it->second << endl;
  }

  return 0;
}
