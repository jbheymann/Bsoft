/**
@file	rwPIR.cpp
@brief	Library routines to read and write PIR sequence files
@author Bernard Heymann
@date	Created: 19990123
@date	Modified: 20250509
**/

#include "rwPIR.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads PIR format sequence files.
@param 	&filename	sequence file name.
@return vector<Bsequence>	set of sequences.

	PIR format:
	>P1;BAC1_HALS1
	BAC1_HALS1, 308 bases, DC74A5E6 checksum.
	 ---------- ---MDPIALT AAVGADLLGD GRPETLWLGI GTLLMLIGTF
	 ASAAD---*

**/
vector<Bsequence>	readPIR(string& filename)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPIR: filename=" << filename << endl;
	
	vector<Bsequence>	seqs;

    // Open seq file read only
	ifstream		fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return seqs;
	}
    
    string			s, seqstr, seqname;
	long			i, nid(0), nseq(0), m(0);
	int				seq_flag(0);
	Bsequence		seq;

//	int 		gaps = 1;								// Retain gaps
//	if ( strstr(molgroup->select, "NOGAP") ) gaps = 0;	// Omit gaps
	
	// Pass to read the data
	nid = nseq = 0;
	while ( !fseq.eof() ) {
		getline(fseq, s);
		if ( s.substr(0,3) == ">P1" ) {
			seqname = s.substr(4);
            seq_flag = 1;
			nid++;
		} else if ( seq_flag == 1 ) {
			seq_flag++;
		} else if ( seq_flag > 1 ) {
			seqstr += s;
            for ( i=0; i<s.length(); i++ ) {
                if ( isalpha(s[i]) ) {
                    seqstr[m] = toupper(s[i]);
                    m++;
                }
                if ( ( s[i] == '-' ) || ( s[i] == '.' ) ) {
//					if ( gaps ) {
	                    seqstr[m] = '-';
    	                m++;
//					}
                }
				if ( s[i] == '*' ) {
					seqstr = seqstr.substr(0,m);
					seq.identifier(seqname);
					seq.sequence(seqstr);
					seqs.push_back(seq);
            		seq_flag = 0;
					seqstr.clear();
					m = 0;
            		nseq++;
				}
            }
		}
    }
    
	fseq.close();
	
	if ( nseq != nid ) {
		cerr << "Number of identifiers (" << nid << ") don't agree with number of sequences (" << nseq << ")" << endl;
		return seqs;
	}
	
    return seqs;
}


/**
@brief 	Writes PIR format sequence files.
@param 	&filename	sequence file name.
@param	seqs	 		set of sequences.
@return int 				number of sequences written (<0 if writing failed).
**/
int			writePIR(string& filename, vector<Bsequence> seqs)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePIR: filename=" << filename << endl;
	
    ofstream        fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}

    int				nseq(0);
    long			j, k, n;
    
 	for ( auto seq: seqs ) {
		n = seq.length();
		fseq << ">P1;" << seq.identifier() << endl;
		fseq << seq.identifier() << ", " << n << " bases, 0 checksum" << endl;
        for ( j=0; j<n; j+=50 ) {
			for ( k=j; k<j+50 && k<n; k+=10 )
				fseq << " " << seq.sequence().substr(k,10);
            fseq << endl;
        }
        fseq << "*" << endl;
        nseq++;
    }
    
    fseq.close();
    
    return nseq;
}
