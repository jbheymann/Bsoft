/**
@file	rwmodel_xyz.h
@brief	Header file for reading and writing Kirkland's xyz atomic model parameters
@author 	Bernard Heymann
@date	Created: 20220426
@date	Modified: 20230706
**/

#include "rwmodel.h"

/* Function prototypes */
Bmodel*		read_model_xyz(vector<string> file_list, map<string,Bcomptype>& atompar);
int			write_model_xyz(string& filename, Bmodel* model, int splt);


