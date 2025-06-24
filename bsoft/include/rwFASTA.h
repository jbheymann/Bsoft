/**
@file	rwFASTA.h
@brief	Header file for reading and writing FASTA sequence files
@author Bernard Heymann
@date	Created: 20001112 
@date	Modified: 20250509

	Format: Protein sequence file format
**/

#include "Bsequence.h"

// I/O prototypes
vector<Bsequence>	readFASTA(string& filename);
int			writeFASTA(string& filename, vector<Bsequence> seqs);
