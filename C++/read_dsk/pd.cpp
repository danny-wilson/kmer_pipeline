// Read a list of aligned, gzipped FASTA DNA/RNA files and output a mean pairwise distance matrix
#include <stdio.h>
#include <map>
#include <iostream>
#include <stdlib.h>
#include <myerror.h>
#include <climits>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <matrix.h>
#include <lotri_matrix.h>
#include <iomanip>

using namespace std;
using namespace myutils;

typedef unsigned long int ULINT;
#define TYPEDEF_ULINT_MAX ULONG_MAX

// Convert ACGTN to 01234, -1 for whitespace, -2 for unrecognized
int encodeACGTN(const char c) {
  switch(c) {
    case 'A':
      return 0;
    case 'C':
      return 1;
    case 'G':
      return 2;
    case 'T': case 'U':
      return 3;
    case 'N': case '?': case '-':
      return 4;
    case ' ': case '\t': case '\r': case '\n': case '\v': case '\b': case '\f':
      return -1;
    default:
      return -2;
  }
}

int main(const int argc, const char* argv[]) {
  if(argc!=3 && argc!=4) {
    cout << "pd is a program for calculating mean pairwise genetic distance for aligned DNA or RNA sequences" << endl;
    cout << "The listed gzipped fasta files are assumed to contain a single label and sequence each" << endl;
    cout << "With two arguments the mean proportion of called sites that differ is output. With three arguments" << endl;
    cout << "the absolute number of called sites that differ is also output." << endl;
    error("Usage: pd fasta.gz_file_list.txt pd_out_filename.txt [abs_pd_out_filename.txt]");
  }
  const char* fagz_listfile = argv[1];
  const char* outfilename = argv[2];
  const char* absfilename = (argc==3) ? "" : argv[3];
  
  // Read the list of gzipped FASTA files
  vector<string> fagz_filename(0);
  ifstream infile(fagz_listfile);
  if(!infile) {
    cout << "Could not open fasta.gz list file " << fagz_listfile << endl;
  }
  int i,j;
  for(i=0;!infile.eof();i++) {
    string fagz_filename_in;
    infile >> fagz_filename_in;
    // Test whether the file exists
    if(fagz_filename_in!="") {
      ifstream fagz_infile(fagz_filename_in.c_str());
      if(!fagz_infile) {
        cout << "Error: Could not open fasta.gz file " << fagz_filename_in << endl;
        exit(13);
      }
      fagz_infile.close();
      fagz_filename.push_back(fagz_filename_in);
    }
  }
  infile.close();
  const int n = fagz_filename.size();
  if(n==0) error("The fasta.gz file list was empty");
  
  // Open the pipes in readiness for reading
  // Vector of pointers to open streams
  vector<FILE*> pin(n,NULL);
  // Maximum number of characters to read at once
  const int pbuff = 10000;
  // Storage for the file buffer
  char buff[pbuff];
  // Indicator of successful read
  char* FGETS;
  for(i=0;i<n;i++) {
    string command = "zcat ";
    command += fagz_filename[i];
    pin[i] = popen(command.c_str(),"r");
    if(pin[i]==NULL) {
      cout << "Error: Could not open stream for file " << i+1 << " " << fagz_filename[i] << endl;
      cout << "Increase the system maximum number of open files if necessary" << endl;
      exit(13);
    }
    // Read the first line, check the formatting, and discard
    FGETS = fgets(buff,pbuff,pin[i]);
    if(FGETS==NULL) {
      cout << "Error: Unexpected end of gzipped fasta file " << fagz_filename[i] << " before label" << endl;
      exit(13);
    }
    // Check the first character is '>'
    if(buff[0]!='>') {
      cout << "Error: Unexpected character \"" << buff[0] << "\" at beginning of gzipped fasta file " << fagz_filename[i] << endl;
      exit(13);
    }
    if(feof(pin[i])) {
      cout << "Error: Unexpectedly reached end of gzipped fasta file " << fagz_filename[i] << "before sequence" << endl;
      exit(13);      
    }
  }
  
  // Create storage for current number of A, C, G, T or N per site
  // Note that ? and - are interpretted as N and U as T
  vector<ULINT> basecount(5);        // Total number of A, C, G, T/U and N/?/-
  Matrix<ULINT> gix(5,n);            // For each base i, gix[i][0..basecount[i]] gives the identity of genomes with that base

  // Create storage for the numerator and denominator of the pairwise difference matrix. Could be very big!
  cout << "Attempting to allocate memory for pairwise distance matrices of size " << n << " by " << n << " genomes" << endl;
  LowerTriangularMatrix<ULINT> num(n,0), den(n,0);
  cout << "Memory successfully allocated" << endl;
  
  // Cycle through sites calculating the number of differences out of mutually called sites
  ULINT pos;
  int gen,nuc;
  for(pos=0;;pos++) {
    if(pos==TYPEDEF_ULINT_MAX) {  // Expected to be at least 4.2 billion unexpected for prokaryotic genomes but may arise for eukaryotic genomes
      cout << "Error: Position " << pos << " exceeds maximum genome size";
      exit(13);
    }
    // Initialize base counts
    for(nuc=0;nuc<5;nuc++) basecount[nuc] = 0;
    // Advance a single base per genome
    for(gen=0;gen<n;gen++) {
      int FGETC = fgetc(pin[gen]);
      // Interpret the base
      int base = encodeACGTN(toupper((char)FGETC));
      while(base==-1) {
        // Repeat until base!=-1 (i.e. skip whitespace)
        FGETC = fgetc(pin[gen]);
        base = encodeACGTN(toupper((char)FGETC));
      }
      // If reached end of file, break
      if(FGETC==EOF) break;
      // If an unrecognized base, throw error
      if(base==-2) {
        cout << "Error: Unrecognized base " << (char)FGETC << " in genome " << fagz_filename[gen] << endl;
        exit(13);
      }
      // Record the index of the current genome in the relevant list
      gix[base][basecount[base]] = gen;
      // Update the number of bases
      basecount[base]++;
    }
    // Test for early exit from loop
    if(gen!=n) {
      // A genome other than the first reached end of file: clearly an error
      if(gen>0) {
        cout << "Error: Genome " << fagz_filename[gen] << " reached end of file at position " << pos << " before genome " << fagz_filename[0] << endl;
        exit(13);
      }
      // First genome reached end of file: test all the others to make sure they reach end of file at the same point
      else {
        for(gen=1;gen<n;gen++) {
          // Advance a single base per genome
          int FGETC = fgetc(pin[gen]);
          // Interpret the base
          int base = encodeACGTN(toupper((char)FGETC));
          while(base==-1) {
            // Repeat until base!=-1 (i.e. skip whitespace)
            FGETC = fgetc(pin[gen]);
            base = encodeACGTN(toupper((char)FGETC));
          }
          // If not reached end of file, error
          if(FGETC!=EOF) {
            cout << "Error: Genome " << fagz_filename[0] << " reached end of file at position " << pos << " before genome " << fagz_filename[gen] << endl;
            exit(13);
          }
        }
        // If reached this point, then cleanly reached end of file simultaneously so quit to final output
        break;
      }
    }
    // Update the pairwise distance calculation
    int nuc2,gen2;
    // If everything is an N, skip what follows
    if(basecount[4]<n) {
      // For invariant sites update the denominator only
      if(basecount[0]==n || basecount[1]==n || basecount[2]==n || basecount[3]==n) {
        for(gen=0;gen<n;gen++) {
          for(gen2=0;gen2<gen;gen2++) {
            //++den.safe(gen,gen2);
			++den[gen][gen2];
          }
        }
      } else {
        // For pairs of genomes involving two non-N characters, increment the denominator and relevant numerator
        for(nuc=0;nuc<4;nuc++) {
          if(basecount[nuc]>0) {
            for(nuc2=0;nuc2<nuc;nuc2++) {
              if(basecount[nuc2]>0) {
                // Update the numerator and denominator for each pair of individuals with these two alleles
                for(gen=0;gen<basecount[nuc];gen++) {
                  const int genome1 = gix[nuc][gen];
                  for(gen2=0;gen2<basecount[nuc2];gen2++) {
                    const int genome2 = gix[nuc2][gen2];
                    //++num.safe(genome1,genome2);
                    //++den.safe(genome1,genome2);
					if(genome1<=genome2) {
					  ++num[genome1][genome2];
					  ++den[genome1][genome2];
					} else {
					  ++num[genome2][genome1];
					  ++den[genome2][genome1];
					}
                  }
                }
              }
            }
            // Update just the denominator for all individuals sharing this allele
            for(gen=0;gen<basecount[nuc];gen++) {
              const int genome1 = gix[nuc][gen];
              for(gen2=0;gen2<gen;gen2++) {
                const int genome2 = gix[nuc][gen2];
                //++den.safe(genome1,genome2);
				if(genome1<=genome2) {
				  ++den[genome1][genome2];
				} else {
				  ++den[genome2][genome1];
				}
              }
            }
          }
        }
        // No need to increment numerator nor denominator for pair involving combinations of Ns and non-Ns or two Ns
      }
    }
    // Loop back round to next position
  }
  
  // Open the output file for writing
  ofstream out(outfilename);
  // Output headers to bip file
  for(i=0;i<n;i++) {
    if(i>0) out << '\t';
    out << fagz_filename[i];
  }
  out << endl;
  // Output mean pairwise distances to file
  out << setprecision(7);
  for(i=0;i<n;i++) {
    for(j=0;j<n;j++) {
      if(j>0) out << '\t';
      if(i==j) out << 0.0;
      else out << double(num.safe(i,j))/double(den.safe(i,j));
    }
    out << endl;
  }
  // Close the output file
  out.close();

  if(absfilename!="") {
    // Open the output file for writing
    ofstream out2(absfilename);
    // Output headers to bip file
    for(i=0;i<n;i++) {
      if(i>0) out2 << '\t';
      out2 << fagz_filename[i];
    }
    out2 << endl;
    // Output absolute pairwise distances to file
    for(i=0;i<n;i++) {
      for(j=0;j<n;j++) {
	if(j>0) out2 << '\t';
	if(i==j) out2 << 0.0;
	else out2 << num.safe(i,j);
      }
      out2 << endl;
    }
    // Close the output file
    out2.close();
  }
  
  // Close the pipes
  for(i=0;i<n;i++) pclose(pin[i]);
  
  return 0;
}
