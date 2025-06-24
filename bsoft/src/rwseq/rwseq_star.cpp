/**
@file	rwseq_star.cpp
@brief	Library routines to read and write sequence files in STAR format
@author 	Bernard Heymann
@date	Created: 19980822 
@date	Modified: 20250509
**/

#include "star.h"
#include "rwseq_star.h"
#include "mol_tags.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads a molecule group from a STAR format file.
@param 	&filename		the file name.
@return vector<Bsequence>	set of sequences.
**/
vector<Bsequence>	read_seq_star(string& filename)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_mol_star: filename=" << filename << endl;

 	Bstar				star;	
	vector<Bsequence>	seqs;

 	if ( star.read(filename) < 0 )
		error_show(filename.c_str(), __FILE__, __LINE__);
	
	if ( star.blocks().size() < 0 ) {
		cerr << "No data blocks found in the STAR file!" << endl;
		return seqs;
	}

	long			nseq(0);
	string			molname;
	Bsequence		seq;

	for ( auto ib: star.blocks() ) {
		molname = ib.tag();
		if ( molname.length() < 1 ) molname = "A";
		seq.identifier(molname);
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG read_mol_star: Molecule " << seq.identifier() << endl;
		seq.sequence(ib.at(MOLECULE_SEQUENCE));
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG read_mol_star: sequence " << seq.sequence() << endl;
		seqs.push_back(seq);
		nseq++;
		molname[0]++;
		if ( molname[0] > 'Z' ) molname[0] = 'A';
	}
		
	if ( nseq < 1 ) nseq--;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_seq_star: " << nseq << " sequences" << endl;
	
	return seqs;
}

/**
@brief 	Writes a molecule group to a STAR format file.
@param 	&filename		the file name.
@param	seqs	 		set of sequences.
@return long 				number of molecules written (<0 if writing failed).
**/
long 		write_seq_star(string& filename, vector<Bsequence> seqs)
{
 	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_seq_star: filename=" << filename << endl;

 	Bstar			star;

//	star.comment();
	star.line_length(120);
	
	long			nseq(0);
	string			chain("A");
	
	for ( auto seq: seqs ) {
		BstarBlock&		block = star.add_block(seq.identifier());
		block[MOLECULE_NAME] = seq.identifier();
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG write_seq_star: sequence " << seq.identifier() << endl;
		block[MOLECULE_LENGTH] = to_string(seq.length());
		block[MOLECULE_SEQUENCE] = seq.sequence();
		nseq++;
		chain[0]++;
	}
			
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_seq_star: " << filename << endl;

	star.write(filename);
	
	return nseq;
}
