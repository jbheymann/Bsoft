/**
@file	rwClustal.cpp
@brief	Library routines to read and write Clustal sequence files
@author Bernard Heymann
@date	Created: 20030309
@date	Modified: 20250509
**/

#include "rwClustal.h"
#include "utilities.h"
#include <fstream>
#include <sstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads Clustal format sequence files.
@param 	&filename		sequence file name.
@return vector<Bsequence>	set of sequences.

	Clustal format:
	
	NP_039940       -----MHAKMNGWAGVRL--VTHCLNTRSRTYVALNMLAFARTPRGVPSCLFNKVWVSRY
	P16720          -----MHAKMNGWAGVRL--VTHCLNTRSRTYVALNMLAFARTPRGVPSCLFNKVWVSRY
	DAA00163        -----MHAKMNGWAGVRL--VTHCLNTRSRTYVALNMLAFARTPRGVPSCLFNKVWVSRY

**/
vector<Bsequence>	readClustal(string& filename)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readClustal: filename=" << filename << endl;
	
	vector<Bsequence>	seqs;

    // Open seq file read only
	ifstream		fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return seqs;
	}
	
	string			s, seqname, seqstr;
    long     		m(0), n(0);
	
	// Read the header line
	getline(fseq, s);
	
	// Read the rest of the data
	while ( !fseq.eof() ) {
		getline(fseq, s);
		if ( s.length() ) {
			istringstream	ss(s);
			ss >> seqname >> seqstr;
			if ( isalnum(seqname[0]) ) {
				if ( m == 0 ) {
					Bsequence 	seq(seqname);
					seq.sequence(seqstr);
					seqs.push_back(seq);
				} else {
					Bsequence&	seq = seqs[n];
					seq.add_sequence(seqstr);
				}
				n++;
			}
			seqname.clear();
		} else {
			if ( n >= seqs.size() ) {
				m++;		// Test for end of first section
				n = 0;
			}
		}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readClustal: sequences read = " << seqs.size() << " lines per sequence = " << m << endl;
		
	fseq.close();
	
    return seqs;
}

/**
@brief 	Writes Clustal format sequence files.
@param 	&filename		sequence file name.
@param	seqs	 		set of sequences.
@return int 				number of sequences written (<0 if writing failed).
**/
int 	writeClustal(string& filename, vector<Bsequence> seqs)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeClustal: filename=" << filename << endl;
	
    ofstream        fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}

    long			m, n, maxlen(0), nseq(0);
    long			max_name_len(10), max_seq_len(60);
	string			molname, s, seq;
    
    for ( auto seq: seqs ) {
		molname = seq.identifier();
		if ( max_name_len < molname.length() ) max_name_len = molname.length();
		n = seq.length();
		if ( maxlen < n ) maxlen = n;
		nseq++;
	}
	if ( max_name_len > 10 ) max_seq_len = 50;
	max_name_len += 6;
	
	time_t		ti = time(NULL);

	fseq << "CLUSTAL alignment file, Written by Bsoft on " << asctime(localtime(&ti)) << endl;

	for ( n = 0; n < maxlen; n += max_seq_len ) {
		fseq << endl;
		for ( auto seq: seqs ) {
			m = max_seq_len;
			if ( maxlen - n < max_seq_len ) m = maxlen - n;
			s.clear();
			if ( n < seq.length() ) s = seq.sequence().substr(n, max_seq_len);
			if ( s.length() < m ) s.append(m-s.length(), '-');
			molname = seq.identifier().substr(0,max_name_len);
			fseq << molname << " " << s << endl;
		}
		fseq << endl;
	}
    
	fseq.close();
    
    return nseq;
}

