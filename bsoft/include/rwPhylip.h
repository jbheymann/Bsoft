/**
@file	rwPhylip.h
@brief	Header file for reading and writing Phylip sequence files
@author Bernard Heymann
@date	Created: 20030308
@date	Modified: 20250509

	Format: Protein sequence file format
**/

#include "Bsequence.h"

// I/O prototypes
vector<Bsequence>	readPhylip(string& filename);
int			writePhylip(string& filename, vector<Bsequence> seqs);
