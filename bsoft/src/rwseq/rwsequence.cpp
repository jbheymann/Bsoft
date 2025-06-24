/**
@file	rwsequence.cpp
@brief	Library routines to read and write sequence files
@author Bernard Heymann
@date	Created: 19980822
@date	Modified: 20250511
**/

#include "seq_util.h"
#include "string_util.h"
#include "utilities.h"

#include "rwseq_star.h"
#include "rwseq_text.h"
#include "rwClustal.h"
#include "rwEMBL.h"
#include "rwFASTA.h"
#include "rwGenBank.h"
#include "rwPhylip.h"
#include "rwPIR.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototypes

/**
@brief 	Reads sequence files.
@param 	&filename			the file name.
@return ivector<Bsequence>		set of sequences.
**/
vector<Bsequence>	read_sequence(string& filename)
{
	string				ext = extension(filename);
	vector<Bsequence>	seqs;

	if ( ext.find("star") != string::npos || ext.find("cif") != string::npos )
    	seqs = read_seq_star(filename);
	else if ( ext.find("txt") != string::npos )
		seqs = read_seq_text(filename);
	else if ( ext.find("aln") != string::npos )
		seqs = readClustal(filename);
	else if ( ext.find("embl") != string::npos )
		seqs = readEMBL(filename);
	else if ( ext.find("fasta") != string::npos )
		seqs = readFASTA(filename);
	else if ( ext.find("gb") != string::npos || ext.find("gp") != string::npos || ext.find("gen") != string::npos )
		seqs = readGenBank(filename);
 	else if ( ext.find("phy") != string::npos )
		seqs = readPhylip(filename);
	else if ( ext.find("pir") != string::npos )
		seqs = readPIR(filename);
	else
		cerr << "Error: File type not supported!" << endl;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Sequences read from file " << filename << ":" << endl;
		for ( auto seq: seqs )
			cout << seq.identifier() << tab << seq.length() << endl;
		cout << endl;
	}
	
	return seqs;
}

/**
@brief 	Reads sequences from a list of files.
@param 	&filename			list of file names.
@return ivector<Bsequence>		set of sequences.
**/
vector<Bsequence>	read_sequence(vector<string>& file_list)
{
	vector<Bsequence>	seqs;
	
	for ( auto f: file_list ) {
		vector<Bsequence>	nseqs = read_sequence(f);
		seqs.insert(seqs.end(), nseqs.begin(), nseqs.end());
	}
	
	return seqs;
}

/**
@brief 	Writes sequence files.
@param 	&filename			the file name.
@param 	seqs				set of sequences.
@return long					number of sequences written, <0 on error.
**/
long		write_sequence(string& filename, vector<Bsequence> seqs)
{
	long				nseq(0);
	string				ext = extension(filename);
	
	if ( verbose & VERB_PROCESS )
		cout << "Writing sequences to " << filename << endl;

	if ( ext.find("star") != string::npos || ext.find("cif") != string::npos )
    	nseq = write_seq_star(filename, seqs);
	else if ( ext.find("txt") != string::npos )
		nseq = write_seq_text(filename, seqs);
	else if ( ext.find("aln") != string::npos )
		nseq = writeClustal(filename, seqs);
	else if ( ext.find("embl") != string::npos )
		nseq = writeEMBL(filename, seqs);
	else if ( ext.find("fasta") != string::npos )
		nseq = writeFASTA(filename, seqs);
	else if ( ext.find("gb") != string::npos || ext.find("gp") != string::npos || ext.find("gen") != string::npos )
		nseq = writeGenBank(filename, seqs);
 	else if ( ext.find("phy") != string::npos )
		nseq = writePhylip(filename, seqs);
	else if ( ext.find("pir") != string::npos )
		nseq = writePIR(filename, seqs);
	else {
		nseq = -1;
		cerr << "Error: File type not supported!" << endl;
	}

	return nseq;
}

long		write_sequence(char *filename, vector<Bsequence> seqs)
{
	string			f(filename);
	
	return write_sequence(f, seqs);
}


