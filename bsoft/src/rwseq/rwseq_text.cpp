/**
@file	rwseq_text.cpp
@brief	Library routines to read and write sequence files in plain text
@author 	Bernard Heymann
@date	Created: 20050419
@date	Modified: 20250509
**/

#include "rwseq_text.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads a molecule group from a text file.
@param 	&filename		the file name.
@return vector<Bsequence>	set of sequences.
**/
vector<Bsequence>	read_seq_text(string& filename)
{
  	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_mol_text: filename=" << filename << endl;
	
	vector<Bsequence>	seqs;

	ifstream    	fseq(filename.c_str());
    if ( fseq.fail() ) {
		cerr << "Error: File " << filename << " not opened!" << endl;
		return seqs;
	}
    
	long			i, nres(0), nna(0), m(0);
	char			na[8] = "ACGTU";
    string			s, seqstr, molname;
	Bsequence		seq;
	
 	while ( !fseq.eof() ) {
		getline(fseq, s);
		if ( s.length() ) {
			if ( isalnum(s[0]) ) {
				seq.identifier(s);
				seqstr.clear();
				m = nna = 0;
			} else {
				s = s.substr(10);
				seqstr += s;
	            for ( i=0; i<s.length(); i++ ) {
    	            if ( isalpha(s[i]) ) {
        	            seqstr[m] = toupper(s[i]);
						if ( strchr(na, seqstr[m]) ) nna++;
						nres++;
                	    m++;
	                }
					if ( ( s[i] == '-' ) || ( s[i] == '.' ) ) {
						seqstr[m] = '-';
						m++;
					}
                }
			}
		} else if ( m ) {
			seqstr = seqstr.substr(0,m);
			if ( 1.5*nna > nres ) seq.type("DNA");
			else seq.type("Protein");
			seqs.push_back(seq);
		}
	}

    fseq.close();	

	return seqs;
}

/**
@brief 	Writes a molecule group to a text file.
@param 	&filename		the file name.
@param	seqs	 		set of sequences.
@return int 				number of molecules written (<0 if writing failed).
**/
int 		write_seq_text(string& filename, vector<Bsequence> seqs)
{
  	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_mol_text: filename=" << filename << endl;
	
    ofstream    	fseq(filename.c_str());

    if ( fseq.fail() ) return -1;

    int     	nseq(0), j, k, n;
	string		seqstr;

	for ( auto seq: seqs ) {
		fseq << seq.identifier() << endl;
		n = seq.length();
        seqstr = seq.sequence();
        for ( j=0; j<n; j+=60 ) {
			fseq << setw(9) << right << j+1;
			for ( k=j; k<j+60 && k<n; k+=10 )
				fseq << " " << seqstr.substr(k,10);
            fseq << endl;
        }
        fseq << endl;
        nseq++;
    }

    fseq.close();
    
	return nseq;
}
