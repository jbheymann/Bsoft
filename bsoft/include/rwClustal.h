/**
@file	rwClustal.h
@brief	Header file for reading and writing Clustal sequence files
@author Bernard Heymann
@date	Created: 20030309
@date	Modified: 20030309

	Format: Protein sequence file format
**/

#include "rwsequence.h"

// I/O prototypes
vector<Bsequence>	readClustal(string& filename);
int 	writeClustal(string& filename, vector<Bsequence> seqs);
