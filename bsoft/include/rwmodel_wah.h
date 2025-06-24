/**
@file	rwmodel_wah.h
@brief	Header file for reading and writing Wayne Hendrickson coordinate files
@author 	Bernard Heymann
@date	Created: 20050217
@date	Modified: 20250318

	Format: Atomic coordinate file format for PROLSQ
**/

#include "rwmodel.h"

// I/O prototypes
Bmodel*	read_model_wah(vector<string> file_list, map<string,Bcomptype>& atompar);
int 	write_model_wah(string& filename, Bmodel* model, int splt);
