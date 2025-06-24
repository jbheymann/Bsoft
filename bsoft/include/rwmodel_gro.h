/**
@file	rwmodel_gro.h
@brief	Header file for reading and writing GROMACS coordinate files
@author 	Bernard Heymann
@date	Created: 19980822
@date	Modified: 20250318

	Format: Atomic coordinate file format for the GROMACS package
**/

#include "rwmodel.h"

// I/O prototypes
Bmodel* 	read_model_gro(vector<string> file_list, map<string,Bcomptype>& atompar);
int 		write_model_gro(string& filename, Bmodel* model, int splt);
