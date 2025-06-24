/**
@file	rwmodel_mol.h
@brief	Header file for reading and writing MDL molfile parameters
@author 	Bernard Heymann
@date	Created: 20060919
@date	Modified: 200250514
**/

#include "rwmodel.h"

/* Function prototypes */
Bmodel*		read_model_mol(vector<string> file_list, map<string,Bcomptype>& atompar);
int			write_model_mol(string& filename, Bmodel* model, int splt);


