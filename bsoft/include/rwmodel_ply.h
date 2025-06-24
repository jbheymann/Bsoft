/**
@file	rwmodel_ply.h
@brief	Header file for reading and writing PLY model parameters
@author 	Bernard Heymann
@date	Created: 20250417
@date	Modified: 20250417
**/

#include "rwmodel.h"

/* Function prototypes */
Bmodel*		read_model_ply(vector<string> file_list);
int			write_model_ply(string& filename, Bmodel* model, int splt);


