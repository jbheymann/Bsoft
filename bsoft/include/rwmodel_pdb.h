/**
@file	rwmodel_pdb.h
@brief	Header file for reading and writing molecular model parameters
@author 	Bernard Heymann
@date	Created: 20211231
@date	Modified: 20230706
**/

#include "rwmodel.h"

/* Function prototypes */
Bmodel*		read_model_pdb(vector<string> file_list);
int			write_model_pdb(string& filename, Bmodel* model, int splt);


