/**
@file	rwFASTA.cpp
@brief	Library routines to read and write FASTA sequence files
@author Bernard Heymann
@date	Created: 19990123
@date	Modified: 20250516
**/

#include "rwFASTA.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads FASTA format sequence files.
@param 	&filename	sequence file name.
@return vector<Bsequence>	set of sequences.

	FASTA format:
	>1hiwa TRIMERIC HIV-1 MATRIX PROTEIN   MOLECULE: HIV-1 MATRIX PROTEIN;   CHAIN:
	VLSGGELDKWEKIRLRPGGKKQYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTIAVLY
	CVHQRIDVKDTKEALDKIEEEQNKSKKKAQQAAAD

**/
vector<Bsequence>	readFASTA(string& filename)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readFASTA: filename=" << filename << endl;
	
	vector<Bsequence>	seqs;

    // Open seq file read only
	ifstream		fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return seqs;
	}
    
    long			i, sl;
	long			m(0), nseq(0);
	string			s, seqstr;
	Bsequence		seq;
	
	while ( !fseq.eof() ) {
		getline(fseq, s);
		sl = s.length();
		if ( sl ) {
			if ( s.back() == '\n' ) {	// Remove trailing newline
				s.pop_back();
				sl = s.length();
			}
			if ( s[0] == '>' ) {
				if ( nseq ) {
					seq.sequence(seqstr.substr(0,m));
					seqs.push_back(seq);
				}
				i = s.find_first_of(" \t\r\n");
				seq = Bsequence(s.substr(1,i-1));
				if ( s.length() > i ) seq.description() = s.substr(i);
				if ( verbose & VERB_FULL )
					cout << seq.identifier() << endl;
				nseq++;
				m = 0;
				seqstr.clear();
			} else {
				seqstr += s;
	            for ( i=0; i<sl; i++ ) {
    	            if ( isalpha(s[i]) ) {
        	            seqstr[m] = toupper(s[i]);
                	    m++;
	                }
					if ( ( s[i] == '-' ) || ( s[i] == '.' ) ) {
						seqstr[m] = '-';
						m++;
					}
                }
            }
		}
    }
    
	if ( seq.length() < 1 ) {		// Clean up last sequence
		seq.sequence(seqstr.substr(0,m));
		seqs.push_back(seq);
	}
    
    fseq.close();
	
    return seqs;
}


/**
@brief 	Writes FASTA format sequence files.
@param 	&filename	sequence file name.
@param	seqs	 		set of sequences.
@return int 				number of sequences written (<0 if writing failed).
**/
int			writeFASTA(string& filename, vector<Bsequence> seqs)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeFASTA: filename=" << filename << endl;
	
    ofstream        fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}

    long			nseq(0), j, n;
	string			seq, s80;
    
	for ( auto seq: seqs ) {
 		n = seq.length();
		fseq << ">" << seq.identifier() << " " << seq.description() << endl;
        for ( j=0; j<n; j+=80 ) {
			s80 = seq.sequence().substr(j,80);
			fseq << s80 << endl;
        }
        nseq++;
    }
    
    fseq.close();
    
    return nseq;
}

