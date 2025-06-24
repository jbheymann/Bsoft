/**
@file	rwPIR.h
@brief	Header file for reading and writing PIR sequence files
@author Bernard Heymann
@date	Created: 19990123
@date	Modified: 20250509
	Format: Protein sequence file format
**/

#include "Bsequence.h"

// I/O prototypes
vector<Bsequence>	readPIR(string& filename);
int			writePIR(string& filename, vector<Bsequence> seqs);
