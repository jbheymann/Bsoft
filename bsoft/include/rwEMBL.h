/**
@file	rwEMBL.h
@brief	Header file for reading and writing EMBL sequence files
@author Bernard Heymann
@date	Created: 19990123
@date	Modified: 20250509

	Format: Protein sequence file format
**/

#include "Bsequence.h"

// I/O prototypes
vector<Bsequence>	readEMBL(string& filename);
int 	writeEMBL(string& filename, vector<Bsequence> seqs);
