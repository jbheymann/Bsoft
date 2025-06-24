/**
@file	rwEMBL.cpp
@brief	Library routines to read and write EMBL sequence files
@author Bernard Heymann
@date	Created: 19990123
@date	Modified: 20250509
**/

#include "rwEMBL.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads EMBL format sequence files.
@param 	&filename		sequence file name.
@return vector<Bsequence>	set of sequences.

	EMBL format:
	ID   AQP1_HUMAN
	DE   AQP1_HUMAN, 527 bases, EB1B8D0 checksum.
	SQ             527 BP
	    ---------- ---------- ---------- ---------- ---------- ----MASEFK

**/
vector<Bsequence>	readEMBL(string& filename)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readEMBL: filename=" << filename << endl;
	
	vector<Bsequence>	seqs;

    // Open seq file read only
	ifstream		fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return seqs;
	}
    
	string			s, tag, seqname, seqstr;
    long			i;
	long			nid(0), nseq(0);
	int				m(0), seq_flag(0);
	Bsequence		seq;

	// Pass to read the data
	while ( !fseq.eof() ) {
		getline(fseq, s);
		tag = s.substr(0,2);
		if ( s.length() > 5 ) s = s.substr(5);
		else s.clear();
		if ( tag == "ID" ) {
			seqname = s;
			nid++;
		} else if ( tag == "//" ) {
			seqstr = seqstr.substr(0,m);
			seq.identifier(seqname);
			seq.sequence(seqstr);
			seqs.push_back(seq);
            seq_flag = 0;
			seqstr.clear();
 			m = 0;
 		} else if ( tag == "SQ" ) {
			seq_flag = 1;
 			nseq++;
        } else if ( seq_flag ) {
			seqstr += s;
			for ( i=0; i<s.length(); ++i ) {
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
    
	fseq.close();

	if ( nid != nseq ) {
		cerr << "Error: Number of identifiers (" << nid << ") don't agree with number of sequences (" << nseq << ")!" << endl;
		return seqs;
	}
	
    return seqs;
}

/**
@brief 	Writes EMBL format sequence files.
@param 	&filename		sequence file name.
@param	seqs	 		set of sequences.
@return int 				number of sequences written (<0 if writing failed).
**/
int 	writeEMBL(string& filename, vector<Bsequence> seqs)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeEMBL: filename=" << filename << endl;
	
    ofstream        fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}

    long			nseq(0), j, k, n;
	string			seqstr;

	for ( auto seq: seqs ) {
		seqstr = seq.sequence();
		n = seqstr.length();
        fseq << "ID   " << seq.identifier() << endl;
        fseq << "DE   " << seq.identifier() << ", " << n << " bases, 0 checksum" << endl;
        fseq << "SQ   " << n << " BP" << endl;
        for ( j=0; j<n; j+=60 ) {
			fseq << "    ";
			for ( k=j; k<j+60 && k<n; k+=10 )
				fseq << " " << seqstr.substr(k,10);
            fseq << endl;
        }
        fseq << "//" << endl;
        nseq++;
    }
    
	fseq.close();
    
    return nseq;
}

