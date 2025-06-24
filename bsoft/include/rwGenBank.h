/**
@file	rwGenBank.h
@brief	Header file for reading and writing GenBank sequence files
@author Bernard Heymann
@date	Created: 20030308
@date	Modified: 20250509

	Format: Protein sequence file format
**/

#include "Bsequence.h"

// I/O prototypes
vector<Bsequence>	readGenBank(string& filename);
int 	writeGenBank(string& filename, vector<Bsequence> seqs);
