/**
@file	rwmol_text.h
@brief	Read and write sequences in plain text
@author 	Bernard Heymann
@date	Created: 20050419
@date	Modified: 20250509
**/

#include "Bsequence.h"

// Function prototypes
vector<Bsequence>	read_seq_text(string& filename);
int 		write_seq_text(string& filename, vector<Bsequence> seqs);

