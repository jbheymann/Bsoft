/**
@file	rwmol_star.h
@brief	Read and write sequences in STAR format
@author 	Bernard Heymann
@date	Created: 19991113
@date	Modified: 20250509
**/

#include "Bsequence.h"

// Function prototypes
vector<Bsequence>	read_seq_star(string& filename);
long 		write_seq_star(string& filename, vector<Bsequence> seqs);

