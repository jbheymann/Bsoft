/**
@file	rwmodel_mol.h
@brief	Header file for reading and writing molecular model parameters
@author Bernard Heymann
@date	Created: 20060919
@date	Modified: 20080408
**/

#include "rwmodel.h"

/* Function prototypes */
Bmodel*		read_model_molecule(vector<string> file_list, string& paramfile);
int			write_model_molecule(string& filename, Bmodel* model);


