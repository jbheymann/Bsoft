/**
@file	rwmodel_cif.h
@brief	Read and write molecules in CIF format
@author 	Bernard Heymann
@date	Created: 19991113
@date	Modified: 20230706
**/

#include "rwmodel.h"

// Function prototypes
Bmodel*		read_model_cif(vector<string> file_list);
long 		write_model_cif(string& filename, Bmodel* model);

