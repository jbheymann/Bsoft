/**
@file	rwPhylip.cpp
@brief	Library routines to read and write Phylip sequence files
@author Bernard Heymann
@date	Created: 20030308
@date	Modified: 20250516
**/

#include "rwPhylip.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads Phylip format sequence files.
@param 	&filename		sequence file name.
@return vector<Bsequence>	set of sequences.

	Phylip format:
	21   1419
	WMBEX6    MT------------------------------------------------
	P89429    MA------------------------------------------------
	P09302    MAEITSLFNNSS--------------------------------------

**/
vector<Bsequence>	readPhylip(string& filename)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPhylip: filename=" << filename << endl;
	
	long				nseq, maxlen;

    // Open seq file read only
	ifstream		fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return vector<Bsequence>();
	}
    	
	fseq >> nseq >> maxlen;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPhylip: nseq=" << nseq << " maxlen=" << maxlen << endl;

	vector<Bsequence>	seqs(nseq);
	string				s, seqname, seqstr;
    long     			m(0), n(0);
	
	// Read the rest of the data
	while ( !fseq.eof() ) {
		getline(fseq, s);
		if ( s.length() ) {
			istringstream	ss(s);
			if ( n < nseq ) {
				ss >> seqname >> seqstr;
				seqs[n].identifier(seqname);
				seqs[n].sequence(seqstr);
				n++;
			} else {
				ss >> seqstr;
				seqs[m].add_sequence(seqstr);
				m++;
				if ( m >= nseq ) m = 0;
			}
		}
	}
	
	fseq.close();
	
    return seqs;
}

/**
@brief 	Writes Phylip format sequence files.
@param 	&filename		sequence file name.
@param	seqs	 		set of sequences.
@return int 				number of sequences written (<0 if writing failed).
**/
int			writePhylip(string& filename, vector<Bsequence> seqs)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePhylip: filename=" << filename << endl;
	
    ofstream        fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}

	bool			interleaved(0);
    long			m, n, maxlen(0);
	int				nseq(seqs.size());
	string			molname, seqstr, s;
 
 	for ( auto seq: seqs )
 		if ( maxlen < seq.length() ) maxlen = seq.length();
   
	fseq << seqs.size() << " " << maxlen << endl;

	if ( interleaved ) {
		for ( n = 0; n < maxlen; n += 50 ) {
			for ( auto seq: seqs ) {
				molname = seq.identifier().substr(0,10);
				seqstr = seq.sequence();
				m = 50;
				if ( maxlen - n < 50 ) m = maxlen - n;
				s.clear();
				if ( n < seq.length() ) s = seqstr.substr(n, 50);
				if ( s.length() < m ) s.append(m-s.length(), '-');
				if ( n == 0 ) fseq << molname << " " << s << endl;
				else fseq << setw(11) << " " << s << endl;
			}
			fseq << endl;
		}
	} else {
		for ( auto seq: seqs ) {
			molname = seq.identifier().substr(0,10);
			seqstr = seq.sequence();
			if ( seqstr.length() < maxlen ) seqstr.append(maxlen-seqstr.length(), '-');
			fseq << molname << " " << seqstr << endl;
		}
		fseq << endl;
	}
    
	fseq.close();
    
    return nseq;
}

